import sys

with open(sys.argv[-1], 'w') as o:
    for filename in sys.argv[1:-1]:
        o.write(filename+'\n')
