## Set os environment
# Author: Yang Yu
# Date: 2018/10/30

import os

# Print the whole environment
print(os.environ)

# Print the path
print(os.environ["PCL_ROOT"])

# Get the value for key PCL_ROOT
os.environ.get("PCL_ROOT")

# Put value in key PCL_ROOT
A = os.environ.get("PCL_ROOT")
os.environ["PCL_ROOT"] = 'abc'
print(os.environ.get("PCL_ROOT"))
os.environ["PCL_ROOT"] = A
print(os.environ.get("PCL_ROOT"))
