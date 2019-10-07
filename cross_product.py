from sympy import *

def Cross(v1, v2):
    if v1.shape == v2.shape:
        if v1.shape != (3,1):
            raise Exception ("shape must be (3,1)")
    else:
        raise Exception ("v1 and v2 must be of the same dimension")
    
    if type(v1) != type(v2):
        raise Exception ("v1 and v2 must be of the same type")
    
    if type(v1) == MutableDenseMatrix or type(v1) == ImmutableDenseMatrix:
        l1 = v1[1]*v2[2] - v1[2]*v2[1]
        l2 = (-1) * (v1[0]*v2[2] - v1[2]*v2[0])
        l3 = v1[0]*v2[1] - v1[1]*v2[0]
        output = Matrix([[l1],[l2],[l3]])
    else:
        l1 = v1[1][0]*v2[2][0] - v1[2][0]*v2[1][0]
        l2 = (-1) * (v1[0][0]*v2[2][0] - v1[2][0]*v2[0][0])
        l3 = v1[0][0]*v2[1][0] - v1[1][0]*v2[0][0]
        output = Matrix([[l1],[l2],[l3]])
    
    return output

if __name__ == "__main__":
    v1 = Matrix([[1],[2],[3]])
    v2 = Matrix([[4],[5],[6]])
    print(Cross(v1,v2))