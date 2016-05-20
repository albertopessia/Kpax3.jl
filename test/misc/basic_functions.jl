# This file is part of Kpax3. License is MIT.

x = [1; 2; 3; 4; 5; 6]
y = [1; 2; 3; 4; 5; 6]

S = 4

N = 1000000

p1 = 1 / factorial(length(x))
p2 = 1 / factorial(S)

v1 = 0
v2 = 0

for i in 1:N
  shuffle!(x)

  if x == [3; 1; 2; 4; 5; 6]
    v1 += 1
  end

  shuffle!(y, S)

  if y == [2; 1; 3; 4; 5; 6]
    v2 += 1
  end
end

v1 /= N
v2 /= N

@test_approx_eq_eps v1 p1 0.005
@test_approx_eq_eps v2 p2 0.005
