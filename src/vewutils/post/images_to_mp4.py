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
    
    # Create a temporary file list for ffmpeg if we have many files
    # For a small number of files, we can use the concat demuxer or pattern
    # For simplicity, we'll use the pattern approach with a temporary directory structure
    # or use the concat demuxer
    
    # Build ffmpeg command using image2 demuxer for each image
    # Use concat filter to combine them into a video
    # Each image is treated as a 1-frame video at the specified framerate
    
    # Create input arguments for each image file
    # Use -framerate to set how long each image should be displayed
    input_args = []
    for img in unique_images:
        abs_path = os.path.abspath(img)
        input_args.extend(['-framerate', str(framerate), '-i', abs_path])
    
    # Build filter_complex to concatenate all inputs and scale
    # Format: [0:v][1:v][2:v]...concat=n=N:v=1,scale=...:a=0[out]
    num_inputs = len(unique_images)
    filter_parts = []
    for i in range(num_inputs):
        filter_parts.append(f'[{i}:v]')
    # Concat filter with scale to ensure even dimensions
    filter_complex = ''.join(filter_parts) + f'concat=n={num_inputs}:v=1,scale=trunc(iw/2)*2:trunc(ih/2)*2[out]'
    
    # Build ffmpeg command
    cmd = [
        'ffmpeg',
    ] + input_args + [
        '-filter_complex', filter_complex,
        '-map', '[out]',
        '-c:v', 'libx264',
        '-pix_fmt', 'yuv420p',
        '-r', str(framerate),  # Output framerate
        '-y',  # Overwrite output file if it exists
        output_file
    ]
    
    print(f"\nRunning ffmpeg command...")
    # Don't print the full command if it's very long (many files)
    if len(unique_images) <= 10:
        print(f"Command: {' '.join(cmd)}")
    else:
        print(f"Command: ffmpeg [with {len(unique_images)} input files] ... {output_file}")
    
    # Run ffmpeg
    result = subprocess.run(
        cmd,
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
