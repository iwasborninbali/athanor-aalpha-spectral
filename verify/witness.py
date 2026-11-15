#!/usr/bin/env python3
"""Exact checker, independent of enumerate.py: SymPy rational linear algebra.
All graphs are declared by explicit edge lists, not decoded from scan masks.
"""
import sympy as S

R=S.Rational

def matrix(n,edges,alpha):
    A=S.zeros(n)
    for u,v in edges: A[u,v]=A[v,u]=1
    d=[sum(A[u,v] for v in range(n)) for u in range(n)]
    M=alpha*S.diag(*d)+(1-alpha)*A
    return A,M,d

def spectra(M):
    return sorted([v for v,m in M.eigenvals().items() for _ in range(m)], key=lambda x:float(x))

def cert(label,assertion,details):
    assert assertion, label
    print('PASS',label,':',details)

triangle=[(0,1),(0,2),(1,2)]
A,M,d=matrix(3,triangle,R(0)); ev=spectra(M)
claimed= -R(2)*2/3
cert('P29 counterexample K3 alpha=0',ev[0]==-1 and -1>claimed, f'lambda_min={ev[0]}, maxcut=2, claimed_upper={claimed}')
cert('P29 corrected Rayleigh formula K3',ev[0]<=R(2)*3/3-R(4)*2/3,'lambda_min=-1 <= corrected_upper=-2/3')

# K1 join K2,2, with a universal vertex 0.
join_edges=[(0,1),(0,2),(0,3),(0,4),(1,3),(1,4),(2,3),(2,4)]
A,M,d=matrix(5,join_edges,R(1,4))
DD=S.diag(*d); r=DD.inv()*M*DD*S.ones(5,1)
cert('P22 informal equality conjecture alpha=1/4',len(set(d))>1 and list(r)==[R(13,4)]*5 and spectra(M)[-1]==R(13,4),f'degrees={d}, D^-1 Aalpha D rowsums={list(r)}, rho=13/4')
A,M,d=matrix(5,join_edges,R(1))
S52_edges=[(u,v) for u in range(5) for v in range(u+1,5) if u<2 or v<2]
_,S52,ds=matrix(5,S52_edges,R(1))
cert('T26(ii) endpoint uniqueness alpha=1',max(d)==max(ds)==4 and sorted(d)!=sorted(ds),f'non-S_5,2 degrees={d}, S_5,2 degrees={ds}, both rho=4; graph is 3-colorable')

# K4 minus edge (2,3): vertices 0 and 1 are interchangeable universal vertices.
minus=[(0,1),(0,2),(0,3),(1,2),(1,3)]
A,M,d=matrix(4,minus,R(1)); ev=spectra(M)
x=S.Matrix([1,0,0,0]); y=S.Matrix([0,1,0,0])
cert('P31(a) endpoint equality',ev[-2]==3 and d!=[3]*4,f'lambda_2={ev[-2]}, degrees={d}, bound=3')
A,M34,d34=matrix(4,minus,R(3,4)); e34=spectra(M34)
cert('P31(a) INTERIOR equality',e34==[2-S.sqrt(2)/2,R(3,2),R(2),2+S.sqrt(2)/2] and e34[-2]==2,
     f'alpha=3/4 spectrum={e34}, lambda_2=2=alpha*n-1, graph != K4')
cert('P16/C17 endpoint automorphism',M*x==3*x and d[0]==d[1] and x[0]!=x[1],f'degrees={d}, x={list(x)} Perron eigenvector and equivalent vertices 0,1 unequal')

# P3: degrees [1,2,1], with no 2-regular component.
path3=[(0,1),(1,2)]
A,M,d=matrix(3,path3,R(0)); adj=spectra(M)
cert('P18 lower equality at alpha=0',adj[-1]==S.sqrt(2),f'lambda(A)={adj[-1]}, connected nonregular degrees={d}')
A,M,d=matrix(3,path3,R(1)); ev=spectra(M)
cert('P18 upper equality at alpha=1',ev[-1]==2 and len(set(d))>1,f'rho={ev[-1]}, Delta=2, connected nonregular')
cert('C33 alpha=1 diameter fails',len(set(ev))==2,f'path diameter 2, only {len(set(ev))} distinct eigenvalues (requires 3)')

A,M,d=matrix(3,triangle,R(1)); ev=spectra(M)
star_bound=(R(1)*(2+1)+S.sqrt((2+1)**2+4*2*(1-2)))/2
cert('P11 endpoint uniqueness alpha=1',ev[-1]==star_bound and len(triangle)==3,f'K3 rho={ev[-1]} equals star_bound={star_bound}, but K3 != K1,2')
A,M,d=matrix(1,[],R(1,4)); ev=spectra(M)
cert('P11/C12 Delta=0 oversight',ev[-1]==0,f'K1 rho=0 but claimed lower bound 1/4')
cert('P8 chromatic r=1 oversight',ev[-1]==0,f'edgeless K1 at alpha=1/4 has PSD zero matrix; alpha < 1/r = 1')

# K3 disjoint K1,3, largest component crossing exactly alpha=1/2.
a=S.symbols('a',real=True)
B=S.Matrix([[3*a,3*(1-a)],[(1-a),a]])
roots=S.solve(B.charpoly().as_expr().subs(B.charpoly().gen,S.Symbol('lambda')) if False else S.det(2*S.eye(2)-B),a)
lam_plus=(4*a+S.sqrt(16*a*a+12*(1-2*a)))/2
cert('Question5 nondifferentiability',roots==[R(1,2)] and lam_plus.subs(a,R(1,2))==2 and S.diff(lam_plus,a).subs(a,R(1,2))==1, f'K3 branch=2, star branch={lam_plus}, crossing at 1/2, left derivative 0, right derivative 1')

# Direct check of paper table: formula for 2x2 bipartite quotient.
def bip_radius(a,b,alpha):
    return (alpha*(a+b)+S.sqrt(alpha**2*(a+b)**2+4*a*b*(1-2*alpha)))/2
for k in range(2,101):
    vals=[bip_radius(3*k+1,3*k+1,R(0)),bip_radius(3,3*k*k,R(0)),bip_radius(1,3*k*k+1,R(0))]
    assert vals[0]==3*k+1 and vals[1]==3*k and vals[2]**2==3*k*k+1
    vals=[bip_radius(3*k+1,3*k+1,R(1,2)),bip_radius(3,3*k*k,R(1,2)),bip_radius(1,3*k*k+1,R(1,2))]
    # Paper actually prints (3k^2+1)/2; exact value is (3k^2+2)/2.
    assert vals==[3*k+1,R(3*k*k+3,2),R(3*k*k+2,2)]
    if k==2:
        cert('p.6 table third entry refuted',vals[2]==7 and R(3*k*k+1,2)==R(13,2),f'k=2 star K1,13: actual lambda(A_1/2)=7; paper table=13/2')
    vals=[bip_radius(3*k+1,3*k+1,R(1)),bip_radius(3,3*k*k,R(1)),bip_radius(1,3*k*k+1,R(1))]
    assert vals==[3*k+1,3*k*k,3*k*k+1]
print('PASS p.6 other table entries and corrected third entry k=2..100, alpha=0,1/2,1 via exact 2x2 quotient')

for a in range(1,11):
    for b in range(1,a+1):
        n=a+b; edges=[(u,v) for u in range(a) for v in range(a,n)]
        for alpha in (R(0),R(1,4),R(1,2),R(3,4),R(1)):
            A,M,d=matrix(n,edges,alpha)
            actual=spectra(M)
            disc=alpha**2*n*n+4*a*b*(1-2*alpha)
            plus=(alpha*n+S.sqrt(disc))/2
            minus=(alpha*n-S.sqrt(disc))/2
            expected=sorted([plus,minus]+[alpha*a]*(b-1)+[alpha*b]*(a-1),key=float)
            assert all(S.simplify(x-y)==0 for x,y in zip(actual,expected)),(a,b,alpha)
print('PASS P37/P38 exact complete bipartite spectra: a=1..10, b=1..a, alpha=0,1/4,1/2,3/4,1; 275 instances')

# Additional checks of the quotient constructions and weight assignments.
for n in range(2,11):
    E=[(u,v) for u in range(n) for v in range(u+1,n)]
    for alpha in (R(0),R(1,4),R(1,2),R(3,4),R(1)):
        A,M,d=matrix(n,E,alpha)
        assert spectra(M)==[alpha*n-1]*(n-1)+[n-1]
print('PASS P36 complete-graph spectra n=2..10 alpha=0,1/4,1/2,3/4,1; 45 instances')

# P14: K2 joined to an independent triple, with independent 2x2 characteristic equation.
join2=[(0,1)]+[(u,v) for u in (0,1) for v in (2,3,4)]
for alpha in (R(0),R(1,4),R(1,2),R(3,4),R(1)):
    A,M,d=matrix(5,join2,alpha)
    Q=S.Matrix([[1+3*alpha,(1-alpha)**2*2*3],[1,0+2*alpha]])
    assert S.simplify(spectra(M)[-1]-spectra(Q)[-1])==0
print('PASS P14 join quotient K2 join E3 at five rational alphas')

for parts in ((2,1,3),(4,2,1,1)):
    n=sum(parts);r=len(parts);offs=[sum(parts[:i]) for i in range(r)]
    edges=[(u,v) for u in range(n) for v in range(u+1,n)
           if not any(o<=u<o+p and o<=v<o+p for o,p in zip(offs,parts))]
    alpha=R(r-1,r)
    A,M,d=matrix(n,edges,alpha)
    x=S.Matrix([R(1,r*p) for p in parts for _ in range(p)])
    assert M*x==R((r-1)*n,r)*x
print('PASS T26(iii) displayed eigenvector weights for partitions (2,1,3) and (4,2,1,1)')