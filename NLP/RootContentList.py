# list the contents of the root directory in HDFS
from snakebite.client import Client

client = Client('localhost', 9000)
for x in client.ls(['/']):
    print x