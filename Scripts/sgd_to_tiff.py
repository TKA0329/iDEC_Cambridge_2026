#!/usr/bin/env python3
"""
Batch-convert Syngene/GeneSnap .sgd gel image files to 16-bit TIFF.

This script extracts the PixelData stream from the OLE/Compound Document
container used by the sample SGD file and writes the original 16-bit grayscale
pixels without rescaling.

Usage:
    python batch_convert_sgd_to_tiff.py /path/to/folder

Optional:
    python batch_convert_sgd_to_tiff.py /path/to/folder --recursive
    python batch_convert_sgd_to_tiff.py /path/to/folder --output /path/to/tiffs

Requires:
    pip install numpy pillow
"""

from __future__ import annotations

import argparse
import struct
import sys
from pathlib import Path

import numpy as np
from PIL import Image

FREESECT = 0xFFFFFFFF
ENDOFCHAIN = 0xFFFFFFFE


def extract_pixeldata_stream(path: Path) -> bytes:
    data = path.read_bytes()

    if data[:8] != bytes.fromhex("D0CF11E0A1B11AE1"):
        raise ValueError("Not an OLE/Compound Document file")

    sector_shift = struct.unpack_from("<H", data, 30)[0]
    sector_size = 1 << sector_shift
    first_dir_sector = struct.unpack_from("<I", data, 48)[0]
    first_difat_sector = struct.unpack_from("<I", data, 68)[0]
    num_difat_sectors = struct.unpack_from("<I", data, 72)[0]
    difat = list(struct.unpack_from("<109I", data, 76))

    def sec(sid: int) -> bytes:
        off = (sid + 1) * sector_size
        return data[off : off + sector_size]

    fat_sids = [x for x in difat if x != FREESECT]
    sid = first_difat_sector
    for _ in range(num_difat_sectors):
        vals = list(struct.unpack("<%dI" % (sector_size // 4), sec(sid)))
        fat_sids.extend(x for x in vals[:-1] if x != FREESECT)
        sid = vals[-1]

    fat = []
    for fsid in fat_sids:
        fat.extend(struct.unpack("<%dI" % (sector_size // 4), sec(fsid)))

    def chain(start: int):
        out = []
        seen = set()
        sid = start
        while (
            sid not in (ENDOFCHAIN, FREESECT)
            and sid < len(fat)
            and sid not in seen
        ):
            out.append(sid)
            seen.add(sid)
            sid = fat[sid]
        return out

    dirdata = b"".join(sec(s) for s in chain(first_dir_sector))

    pixel_start = None
    pixel_size = None
    for i in range(0, len(dirdata), 128):
        entry = dirdata[i : i + 128]
        if len(entry) < 128:
            break
        name_len = struct.unpack_from("<H", entry, 64)[0]
        if name_len >= 2:
            name = entry[: name_len - 2].decode("utf-16le", "replace")
        else:
            name = ""

        if name == "PixelData":
            pixel_start = struct.unpack_from("<I", entry, 116)[0]
            pixel_size = struct.unpack_from("<Q", entry, 120)[0]
            break

    if pixel_start is None or pixel_size is None:
        raise ValueError("PixelData stream not found")

    return b"".join(sec(s) for s in chain(pixel_start))[:pixel_size]


def convert_one(src: Path, dst: Path) -> tuple[int, int, int, int]:
    pixel = extract_pixeldata_stream(src)

    if len(pixel) < 26:
        raise ValueError("PixelData stream is too short")

    image_type, width, height, a, b, c, d = struct.unpack_from("<H6I", pixel, 0)
    expected = width * height * 2

    if width <= 0 or height <= 0:
        raise ValueError(f"Invalid image dimensions: {width} x {height}")

    if len(pixel) < 26 + expected:
        raise ValueError(
            f"PixelData shorter than expected for {width} x {height} 16-bit image"
        )

    pixels = np.frombuffer(
        pixel, dtype="<u2", offset=26, count=width * height
    ).reshape(height, width)

    dst.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(pixels).save(dst, compression="raw")

    return width, height, int(pixels.min()), int(pixels.max())


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Batch-convert Syngene/GeneSnap .sgd files to 16-bit TIFF"
    )
    parser.add_argument("input", type=Path, help="Folder containing .sgd files")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output folder. Default: a 'converted_tiff' folder inside the input folder",
    )
    parser.add_argument(
        "--recursive", action="store_true", help="Search subfolders recursively"
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing TIFF files"
    )
    args = parser.parse_args()

    input_dir = args.input.expanduser().resolve()
    if not input_dir.is_dir():
        print(f"ERROR: Not a folder: {input_dir}", file=sys.stderr)
        return 2

    output_dir = (
        args.output.expanduser().resolve()
        if args.output
        else input_dir / "converted_tiff"
    )

    pattern = "**/*.sgd" if args.recursive else "*.sgd"
    files = sorted(input_dir.glob(pattern))

    if not files:
        print(f"No .sgd files found in {input_dir}")
        return 0

    print(f"Found {len(files)} .sgd file(s)")
    print(f"Output folder: {output_dir}")
    print()

    ok = 0
    failed = 0

    for src in files:
        if args.recursive:
            rel = src.relative_to(input_dir).with_suffix(".tif")
            dst = output_dir / rel
        else:
            dst = output_dir / f"{src.stem}.tif"

        if dst.exists() and not args.overwrite:
            print(f"SKIP  {src.name} -> {dst.name} (already exists)")
            continue

        try:
            width, height, vmin, vmax = convert_one(src, dst)
            print(
                f"OK    {src.name} -> {dst.name}  "
                f"[{width}x{height}, 16-bit, range {vmin}-{vmax}]"
            )
            ok += 1
        except Exception as e:
            print(f"FAIL  {src.name}: {e}")
            failed += 1

    print()
    print(f"Finished: {ok} converted, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())