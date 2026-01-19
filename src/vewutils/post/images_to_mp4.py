#!/usr/bin/env python3
"""
Module for converting image files to MP4 video using ffmpeg.
"""

import argparse
import subprocess
import sys
import os
import glob
from pathlib import Path
from typing import List


def images_to_mp4(
    image_files: List[str],
    output_file: str,
    framerate: float = 30.0
) -> None:
    """
    Convert image files to MP4 video using ffmpeg.
    
    Args:
        image_files: List of image file paths or glob patterns
        output_file: Path to output MP4 file
        framerate: Frame rate for the output video (fps)
    """
    # Expand glob patterns and collect all image files
    all_images = []
    for pattern in image_files:
        # Check if it's a glob pattern
        if '*' in pattern or '?' in pattern or '[' in pattern:
            matched = sorted(glob.glob(pattern))
            if not matched:
                print(f"Warning: No files matched pattern: {pattern}")
            all_images.extend(matched)
        else:
            # Regular file path
            if os.path.exists(pattern):
                all_images.append(pattern)
            else:
                print(f"Warning: File not found: {pattern}")
    
    if not all_images:
        raise ValueError("No image files found. Please check your input files/patterns.")
    
    # Remove duplicates while preserving order
    seen = set()
    unique_images = []
    for img in all_images:
        if img not in seen:
            seen.add(img)
            unique_images.append(img)
    
    print(f"Found {len(unique_images)} image file(s)")
    print(f"Output file: {output_file}")
    print(f"Frame rate: {framerate} fps")
    
    # Determine image extension from first file
    first_ext = os.path.splitext(unique_images[0])[1] or '.png'
    
    # Create temporary directory with numbered symlinks for ffmpeg pattern matching
    import tempfile
    import shutil
    temp_dir = None
    pattern_path = None
    
    try:
        # Check if all files are in the same directory and follow a simple pattern
        # If so, we might be able to use glob directly
        all_same_dir = len(set(os.path.dirname(os.path.abspath(img)) for img in unique_images)) == 1
        
        if all_same_dir and len(image_files) == 1 and ('*' in image_files[0] or '?' in image_files[0]):
            # Can use glob pattern directly
            work_dir = os.path.dirname(os.path.abspath(unique_images[0]))
            pattern = os.path.basename(image_files[0])
            pattern_path = os.path.join(work_dir, pattern)
        else:
            # Create temporary directory with numbered symlinks
            temp_dir = tempfile.mkdtemp(prefix='ffmpeg_images_')
            
            # Create numbered symlinks: img_0000.png, img_0001.png, etc.
            num_digits = len(str(len(unique_images) - 1))
            pattern = f'img_%0{num_digits}d{first_ext}'
            
            for idx, img in enumerate(unique_images):
                # Format with proper zero-padding
                link_name = f'img_{idx:0{num_digits}d}{first_ext}'
                link_path = os.path.join(temp_dir, link_name)
                os.symlink(os.path.abspath(img), link_path)
            
            # For pattern_type sequence (not glob), use the pattern directly
            pattern_path = pattern
        
        # Build ffmpeg command using pattern-based approach
        # Use pad instead of scale to maintain original dimensions while ensuring even size
        if temp_dir:
            # Use sequence pattern for numbered symlinks
            pattern_type = 'sequence'
            input_pattern = pattern_path
            work_dir = temp_dir
        else:
            # Use glob pattern for direct file matching
            pattern_type = 'glob'
            input_pattern = pattern_path
            work_dir = os.path.dirname(pattern_path) if os.path.dirname(pattern_path) else '.'
        
        # Use absolute path for output file
        abs_output = os.path.abspath(output_file)
        
        cmd = [
            'ffmpeg',
            '-framerate', str(framerate),
            '-pattern_type', pattern_type,
            '-i', input_pattern,
            '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2,format=yuv420p',
            '-c:v', 'libx264',
            '-pix_fmt', 'yuv420p',
            '-crf', '23',
            '-movflags', '+faststart',
            '-y',  # Overwrite output file if it exists
            abs_output
        ]
        
        print(f"\nRunning ffmpeg command...")
        print(f"Command: {' '.join(cmd)}")
        
        # Run ffmpeg in the working directory
        result = subprocess.run(
            cmd,
            cwd=work_dir,
            capture_output=True,
            text=True,
            check=False
        )
        
        if result.returncode != 0:
            print(f"Error: ffmpeg failed with return code {result.returncode}")
            print(f"stderr: {result.stderr}")
            raise RuntimeError(f"ffmpeg command failed: {result.stderr}")
        
        print(f"\nSuccessfully created video: {output_file}")
        if result.stdout:
            print(f"ffmpeg output: {result.stdout}")
    
    finally:
        # Clean up temporary directory
        if temp_dir and os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)


def get_parser():
    parser = argparse.ArgumentParser(
        add_help=False,
        description="Convert image files to MP4 video using ffmpeg"
    )
    parser.add_argument(
        "image_files",
        nargs='+',
        help="Image file paths or glob patterns (e.g., '*.png', 'frame_*.jpg')"
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        dest="output",
        help="Path to output MP4 file"
    )
    parser.add_argument(
        "--framerate",
        type=float,
        default=30.0,
        help="Frame rate for the output video in fps (default: 30.0)"
    )
    return parser


def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    
    print(f"Converting images to MP4 video...")
    print(f"Input files/patterns: {args.image_files}")
    
    try:
        images_to_mp4(
            args.image_files,
            args.output,
            framerate=args.framerate
        )
        print("Conversion completed successfully.")
        return 0
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
