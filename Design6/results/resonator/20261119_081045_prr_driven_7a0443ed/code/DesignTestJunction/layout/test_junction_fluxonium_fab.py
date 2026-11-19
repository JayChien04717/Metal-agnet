"""Generate the vertical F test coupon, fabrication GDS and reference."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build import main


if __name__ == '__main__':
    main('fluxonium')
