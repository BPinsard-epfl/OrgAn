import os
import sys
current_dir = os.getcwd()
target_dir_relative = os.path.join("src", "OrgAn")
target_dir_absolute = os.path.abspath(os.path.join(current_dir, target_dir_relative))
sys.path.append(target_dir_absolute)

from chromato import get_elution_order
from functions import gives_data_frame

solutes = gives_data_frame("tests/test_data.csv")
print(get_elution_order(solutes))