print "ply"
print "format ascii 1.0"
numFaces = int(raw_input())
numVertex = 0
numTri = 0
n, m, P = [], [], []
for k in range(numFaces):
    _n, _m = map(int, raw_input().split())
    p = []
    for i in range(_n + 1):
        for j in range(_m + 1):
            _p = map(float, raw_input().split())
            p.append(_p)
    n.append(_n)
    m.append(_m)
    numVertex += (_n + 1) * (_m + 1)
    numTri += _n * _m * 2
    P.append(p)
print "element vertex %d" % numVertex
print "property float x"
print "property float y"
print "property float z"
print "element face %d" % numTri
print "property list uchar int vertex_indices"
print "end_header"
for p in P:
    for _p in p:
        print "%.6lf %.6lf %.6lf" % (_p[0], _p[1], _p[2])
cur = 0
for k in range(numFaces):
    for i in range(_n):
        for j in range(_m):
            print "3 %d %d %d" % (
                cur + i * (_m + 1) + j,
                cur + (i + 1) * (_m + 1) + j,
                cur + i * (_m + 1) + j + 1
            )
            print "3 %d %d %d" % (
                cur + (i + 1) * (_m + 1) + (j + 1),
                cur + (i + 1) * (_m + 1) + j,
                cur + i * (_m + 1) + j + 1
            )
    cur += (_n + 1) * (_m + 1)
