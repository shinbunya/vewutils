"""Scan for newly-completed ADCIRC forecast cycles and generate/upload reports.

Meant to be run periodically by cron (e.g. every 10 minutes); this script
does a single scan-and-process pass and exits, rather than running as a
daemon. See [watch] and [sftp] in the report TOML config for settings.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import selectors
import subprocess
import sys
import time
from collections import deque
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

DEFAULT_MAX_AGE_DAYS = 3
DEFAULT_QUIET_SECONDS = 60
DEFAULT_MAX_ATTEMPTS = 5
DEFAULT_BACKOFF_MINUTES = [0, 30, 60, 120, 240]
STATE_FILENAME = '.generate_report_state.json'
LOG_FILENAME = '.generate_report_log.txt'
REPORT_SUBDIR = 'report'
GENERATE_REPORT_TIMEOUT_SECONDS = 60 * 60

# Matches the cycle date/hour embedded in an archive cycle path, e.g.
# ".../archive/20260723/hour_12/adcirc/forecast/forecast_base".
ARCHIVE_PATH_RE = re.compile(r'/(\d{8})/hour_(\d{2})(?:/|$)')


def parse_cycle_datetime(path: Path) -> datetime | None:
    """Parse a cycle directory's YYYYMMDD/hour_HH into a datetime, or None."""
    match = ARCHIVE_PATH_RE.search(str(path) + '/')
    if not match:
        return None
    ymd, hh = match.groups()
    try:
        return datetime.strptime(f'{ymd}{hh}', '%Y%m%d%H')
    except ValueError:
        return None


def newest_mtime(path: Path) -> float:
    """Most recent mtime (epoch seconds) among all files under path, recursive."""
    newest = path.stat().st_mtime
    for root, _dirs, files in os.walk(path):
        for name in files:
            try:
                mtime = os.path.getmtime(os.path.join(root, name))
            except OSError:
                continue
            newest = max(newest, mtime)
    return newest


def _load_state(state_path: Path) -> dict[str, Any] | None:
    if not state_path.is_file():
        return None
    try:
        with open(state_path, encoding='utf-8') as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        return None


def _save_state(state_path: Path, attempts: int, error: str) -> None:
    with open(state_path, 'w', encoding='utf-8') as f:
        json.dump({
            'attempts': attempts,
            'last_attempt': datetime.utcnow().isoformat() + 'Z',
            'last_error': error.strip()[-2000:],
        }, f, indent=2)


def _backoff_ready(state: dict[str, Any], backoff_minutes: list[int]) -> bool:
    """True once enough time has passed since last_attempt for another try."""
    attempts = state.get('attempts', 0)
    if attempts <= 0:
        return True
    wait_minutes = backoff_minutes[min(attempts - 1, len(backoff_minutes) - 1)]
    try:
        last_attempt = datetime.fromisoformat(str(state['last_attempt']).rstrip('Z'))
    except (KeyError, ValueError):
        return True
    return datetime.utcnow() >= last_attempt + timedelta(minutes=wait_minutes)


def resolve_cutoff(
        *,
        min_cycle_date: str | None,
        min_cycle_hour: int | None,
        max_age_days: float) -> datetime:
    """Resolve the cutoff below which cycles are ignored.

    If min_cycle_date is given, the cutoff is that absolute date (+hour,
    default 0) -- a fixed watermark that doesn't move as time passes. Otherwise
    falls back to the relative `now - max_age_days` window.
    """
    if min_cycle_date:
        try:
            cutoff_date = datetime.strptime(min_cycle_date, '%Y-%m-%d')
        except ValueError as exc:
            raise ValueError(
                f'min_cycle_date must be YYYY-MM-DD, got {min_cycle_date!r}'
            ) from exc
        return cutoff_date + timedelta(hours=int(min_cycle_hour or 0))
    return datetime.utcnow() - timedelta(days=max_age_days)


def select_candidates(pattern: str, cutoff: datetime) -> list[Path]:
    """Glob pattern, keep directories whose parsed cycle time is >= cutoff."""
    dated: list[tuple[datetime, Path]] = []
    for match in glob.glob(pattern):
        path = Path(match)
        if not path.is_dir():
            continue
        cycle_dt = parse_cycle_datetime(path)
        if cycle_dt is None or cycle_dt < cutoff:
            continue
        dated.append((cycle_dt, path))
    dated.sort(key=lambda item: item[0])
    return [path for _, path in dated]


def is_ready(
        candidate: Path,
        *,
        quiet_seconds: float,
        max_attempts: int,
        backoff_minutes: list[int]) -> tuple[bool, str]:
    """Return (ready, reason); reason explains why not, when ready is False."""
    report_dir = candidate / REPORT_SUBDIR
    state_path = candidate / STATE_FILENAME
    state = _load_state(state_path)

    if report_dir.is_dir() and state is None:
        return False, 'already has report/ and no pending retry state'

    if state is not None:
        attempts = state.get('attempts', 0)
        if attempts >= max_attempts:
            return False, f'gave up after {attempts} failed attempt(s)'
        if not _backoff_ready(state, backoff_minutes):
            return False, f'backing off ({attempts} attempt(s) so far)'

    age_seconds = time.time() - newest_mtime(candidate)
    if age_seconds < quiet_seconds:
        return False, f'still being written (newest file {age_seconds:.0f}s old)'

    return True, ''


LOG_TAIL_LINES = 200


def run_generate_report(config_path: str, candidate: Path) -> tuple[bool, str]:
    """Run generate-report for one candidate directory as a subprocess.

    generate-report's own stdout/stderr (merged, since separating them would
    lose chronological ordering) are streamed live to
    <candidate>/.generate_report_log.txt as they're produced, line by line --
    rather than captured in full and only written after the subprocess exits
    -- so `tail -f` on that file shows real-time progress for a report that's
    still generating, not just the previous attempt's full output once this
    one finishes.

    Deliberately kept outside report/ (a sibling of .generate_report_state.json
    in the candidate directory itself) so it never gets swept into the SFTP
    upload, which only publishes report/'s contents.
    """
    cmd = [
        sys.executable, '-m', 'vewutils.cli', 'plot', 'generate-report',
        '--config', str(config_path),
        '--data-dir', str(candidate),
        '--skip-existing', '--skip-on-error',
    ]
    log_path = candidate / LOG_FILENAME
    timestamp = datetime.utcnow().isoformat() + 'Z'
    tail: deque[str] = deque(maxlen=LOG_TAIL_LINES)

    with open(log_path, 'a', encoding='utf-8') as log_file:
        log_file.write(f'\n=== {timestamp} starting: {" ".join(cmd)} ===\n')
        log_file.flush()

        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1,
        )
        sel = selectors.DefaultSelector()
        sel.register(proc.stdout, selectors.EVENT_READ)
        deadline = time.monotonic() + GENERATE_REPORT_TIMEOUT_SECONDS
        timed_out = False
        try:
            while True:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    timed_out = True
                    break
                if sel.select(timeout=min(remaining, 1.0)):
                    line = proc.stdout.readline()
                    if line == '' and proc.poll() is not None:
                        break
                    if line:
                        log_file.write(line)
                        log_file.flush()
                        tail.append(line)
                elif proc.poll() is not None:
                    break
        finally:
            sel.close()

        if timed_out:
            proc.kill()
            proc.wait()
            log_file.write(f'=== {timestamp} timed out after {GENERATE_REPORT_TIMEOUT_SECONDS}s ===\n')
            return False, f'timed out after {GENERATE_REPORT_TIMEOUT_SECONDS}s'

        returncode = proc.wait()
        log_file.write(f'=== {timestamp} exit code {returncode} ===\n')

    if returncode != 0:
        return False, ''.join(tail) or f'exit code {returncode}'
    return True, ''


def render_remote_path(template: str, cycle_dt: datetime) -> str:
    """Substitute YYYYMMDDHH before YYYY, since YYYY is a prefix of YYYYMMDDHH."""
    rendered = template.replace('YYYYMMDDHH', cycle_dt.strftime('%Y%m%d%H'))
    rendered = rendered.replace('YYYY', cycle_dt.strftime('%Y'))
    return rendered


def _sftp_makedirs(sftp, remote_dir: str) -> None:
    """mkdir -p equivalent over SFTP."""
    current = ''
    for part in remote_dir.strip('/').split('/'):
        current += '/' + part
        try:
            sftp.stat(current)
        except FileNotFoundError:
            sftp.mkdir(current)


def _sftp_upload_dir(sftp, local_dir: Path, remote_dir: str) -> None:
    """Recursively upload local_dir's contents into remote_dir over SFTP."""
    _sftp_makedirs(sftp, remote_dir)
    for item in sorted(local_dir.iterdir()):
        remote_item = f'{remote_dir.rstrip("/")}/{item.name}'
        if item.is_dir():
            _sftp_upload_dir(sftp, item, remote_item)
        else:
            sftp.put(str(item), remote_item)


def upload_report(
        sftp_cfg: dict[str, Any],
        report_dir: Path,
        cycle_dt: datetime) -> tuple[bool, str]:
    """Upload report_dir to the remote host/path described by sftp_cfg."""
    import paramiko

    host = sftp_cfg.get('host')
    username = sftp_cfg.get('username')
    key_path = sftp_cfg.get('key_path')
    remote_path_template = sftp_cfg.get('remote_path_template')
    if not (host and username and key_path and remote_path_template):
        return False, (
            'sftp config incomplete: host, username, key_path, and '
            'remote_path_template are all required'
        )
    port = int(sftp_cfg.get('port', 22))
    key_path = os.path.expanduser(str(key_path))

    remote_dir = render_remote_path(remote_path_template, cycle_dt)

    client = paramiko.SSHClient()
    client.load_system_host_keys()
    try:
        client.connect(host, port=port, username=username, key_filename=key_path)
        sftp = client.open_sftp()
        try:
            _sftp_upload_dir(sftp, report_dir, remote_dir)
        finally:
            sftp.close()
    except Exception as exc:
        return False, f'sftp upload to {host}:{remote_dir} failed: {exc}'
    finally:
        client.close()

    print(f'Uploaded {report_dir} to {host}:{remote_dir}')
    return True, ''


def scan_and_process(
        config_path: str,
        *,
        pattern: str | None = None,
        max_age_days: float | None = None,
        min_cycle_date: str | None = None,
        min_cycle_hour: int | None = None,
        quiet_seconds: float | None = None,
        max_attempts: int | None = None) -> int:
    """Single scan-and-process pass over all cycles matching [watch].pattern."""
    import tempfile

    from vewutils.plot.generate_report import load_report_config

    # A watch template's own [report].data_dir is just a placeholder -- the
    # real value is supplied per-candidate via run_generate_report's
    # --data-dir. Override it here with something that always exists so
    # load_report_config's validation doesn't reject the template for lacking
    # a real data_dir; only [watch]/[sftp] are actually read from this call.
    config = load_report_config(config_path, data_dir_override=tempfile.gettempdir())
    watch_cfg = config['watch']
    sftp_cfg = config['sftp']

    effective_pattern = pattern or watch_cfg.get('pattern')
    if not effective_pattern:
        print(
            'Error: a scan pattern is required (set [watch].pattern in the '
            'config, or pass --pattern)',
            file=sys.stderr,
        )
        return 1

    effective_max_age_days = (
        max_age_days if max_age_days is not None
        else watch_cfg.get('max_age_days', DEFAULT_MAX_AGE_DAYS)
    )
    effective_min_cycle_date = (
        min_cycle_date if min_cycle_date is not None
        else watch_cfg.get('min_cycle_date')
    )
    effective_min_cycle_hour = (
        min_cycle_hour if min_cycle_hour is not None
        else watch_cfg.get('min_cycle_hour', 0)
    )
    effective_quiet_seconds = (
        quiet_seconds if quiet_seconds is not None
        else watch_cfg.get('quiet_seconds', DEFAULT_QUIET_SECONDS)
    )
    effective_max_attempts = (
        max_attempts if max_attempts is not None
        else watch_cfg.get('max_attempts', DEFAULT_MAX_ATTEMPTS)
    )
    backoff_minutes = watch_cfg.get('backoff_minutes', DEFAULT_BACKOFF_MINUTES)

    cutoff = resolve_cutoff(
        min_cycle_date=effective_min_cycle_date,
        min_cycle_hour=effective_min_cycle_hour,
        max_age_days=effective_max_age_days,
    )
    candidates = select_candidates(effective_pattern, cutoff)
    print(f'Found {len(candidates)} candidate cycle(s) at or after {cutoff}')

    for candidate in candidates:
        ready, reason = is_ready(
            candidate,
            quiet_seconds=effective_quiet_seconds,
            max_attempts=effective_max_attempts,
            backoff_minutes=backoff_minutes,
        )
        if not ready:
            print(f'Skip {candidate}: {reason}')
            continue

        state_path = candidate / STATE_FILENAME
        print(f'Processing {candidate}')
        success, error = run_generate_report(config_path, candidate)
        if not success:
            attempts = (_load_state(state_path) or {}).get('attempts', 0) + 1
            _save_state(state_path, attempts, error)
            last_line = error.strip().splitlines()[-1] if error.strip() else error
            print(f'  generate-report failed (attempt {attempts}): {last_line}')
            continue

        if sftp_cfg:
            cycle_dt = parse_cycle_datetime(candidate)
            uploaded, upload_error = upload_report(
                sftp_cfg, candidate / REPORT_SUBDIR, cycle_dt
            )
            if not uploaded:
                attempts = (_load_state(state_path) or {}).get('attempts', 0) + 1
                _save_state(state_path, attempts, upload_error)
                print(f'  upload failed (attempt {attempts}): {upload_error}')
                continue

        if state_path.is_file():
            state_path.unlink()
        print('  done')

    return 0


def get_parser():
    parser = argparse.ArgumentParser(
        add_help=False,
        description=(
            'Scan for newly-completed ADCIRC forecast cycles and run '
            'generate-report (optionally uploading the result via SFTP) on '
            'each one. Meant to be invoked periodically by cron.'
        ),
    )
    parser.add_argument(
        '--config',
        required=True,
        help='Path to report TOML config (also provides [watch]/[sftp] settings)',
    )
    parser.add_argument(
        '--pattern',
        default=None,
        help='Glob pattern for cycle directories (overrides [watch].pattern)',
    )
    parser.add_argument(
        '--max-age-days',
        type=float,
        default=None,
        help=(
            'Ignore cycles older than this many days from now (overrides '
            f'[watch].max_age_days, default {DEFAULT_MAX_AGE_DAYS}). Ignored '
            'if --min-cycle-date/[watch].min_cycle_date is set.'
        ),
    )
    parser.add_argument(
        '--min-cycle-date',
        default=None,
        metavar='YYYY-MM-DD',
        help=(
            'Ignore cycles before this absolute date (overrides '
            '[watch].min_cycle_date). A fixed watermark instead of a rolling '
            '--max-age-days window; takes precedence over --max-age-days when set.'
        ),
    )
    parser.add_argument(
        '--min-cycle-hour',
        type=int,
        default=None,
        metavar='HH',
        help=(
            'Hour (0-23) combined with --min-cycle-date (overrides '
            '[watch].min_cycle_hour, default 0)'
        ),
    )
    parser.add_argument(
        '--quiet-seconds',
        type=float,
        default=None,
        help=(
            'Required idle time since a cycle\'s newest file before '
            f'processing it (overrides [watch].quiet_seconds, default {DEFAULT_QUIET_SECONDS})'
        ),
    )
    parser.add_argument(
        '--max-attempts',
        type=int,
        default=None,
        help=f'Give up on a cycle after this many failed attempts (overrides [watch].max_attempts, default {DEFAULT_MAX_ATTEMPTS})',
    )
    return parser


def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    return scan_and_process(
        args.config,
        pattern=args.pattern,
        max_age_days=args.max_age_days,
        min_cycle_date=args.min_cycle_date,
        min_cycle_hour=args.min_cycle_hour,
        quiet_seconds=args.quiet_seconds,
        max_attempts=args.max_attempts,
    )


if __name__ == '__main__':
    sys.exit(main())
