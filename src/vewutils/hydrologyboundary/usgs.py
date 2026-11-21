def getQH(siteids,staDT,endDT,cachedir=None):
    """Return USGS discharge (m^3/s) and gage height (m) for each site.

    Uses the USGS Water Data API (continuous values). staDT and endDT
    should be timezone-aware; naive times are taken as UTC. Returns four
    lists, one entry per site: discharge times, discharge values, gage
    height times and gage height values. Times are timezone-aware (UTC).
    """
    from datetime import datetime
    import hashlib
    import os
    from pathlib import Path
    import json

    from vewutils.utils import usgs_waterdata

    if type(siteids) != list:
        siteids = [siteids]

    stStaDT = staDT.strftime("%Y-%m-%dT%H:%M%z")
    stEndDT = endDT.strftime("%Y-%m-%dT%H:%M%z")

    timqs = []
    valqs = []
    timhs = []
    valhs = []

    for siteid in siteids:
        cache_key = 'waterdata|{}|00060,00065|{}|{}'.format(siteid, stStaDT, stEndDT)
        if cachedir:
            hasher = hashlib.md5()
            hasher.update(cache_key.encode('utf-8'))
            cachefile = os.path.join(cachedir, hasher.hexdigest()+'.json')
        else:
            cachefile = None
        print('Target = USGS {} 00060,00065 {} to {}'.format(siteid, stStaDT, stEndDT))
        if cachefile and Path(cachefile).is_file():
            print('Found in local cache: {:s}. Loading.'.format(cachefile))
            with open(cachefile, 'r') as f:
                json_res = json.load(f)
        else:
            print('Not found in local cache. Downloading.')
            series = usgs_waterdata.get_continuous(siteid, ['00060', '00065'], staDT, endDT)
            for code in ('00060', '00065'):
                if code not in series:
                    raise RuntimeError(
                        'No USGS {} data for site {} from {} to {}'.format(code, siteid, stStaDT, stEndDT)
                    )
            json_res = {
                code: {
                    'time': [t.isoformat() for t in series[code].index],
                    'value': series[code].tolist(),
                }
                for code in ('00060', '00065')
            }

            if cachefile:
                with open(cachefile, 'w') as f:
                    json.dump(json_res, f)

        # Discharge
        timq = [datetime.fromisoformat(t) for t in json_res['00060']['time']]
        valq = [v*0.3048**3 for v in json_res['00060']['value']]

        # Gage height
        timh = [datetime.fromisoformat(t) for t in json_res['00065']['time']]
        valh = [v*0.3048 for v in json_res['00065']['value']]

        # Store
        timqs.append(timq)
        valqs.append(valq)
        timhs.append(timh)
        valhs.append(valh)

    return timqs, valqs, timhs, valhs
