import os
import sys, pathlib

pth = pathlib.Path(__file__).parent
sys.path.insert(0, str(pth))
os.environ["DSR_DIR"] = str(pth)
