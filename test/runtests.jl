using Echidna
using Base.Test

# write your own tests here
lb = 4.0
ub = 9.0
mr = MOGA_Real(lb, ub)
res = [rand(mr) for i in 1:1000]
@test all(r>lb for r in res) && all(r<ub for r in res)
