# This file is part of Kpax3. License is MIT.

support = KSupport(10, 6, 1, 1);

R1 = [1; 1; 2; 2; 2; 3]
R2 = [1; 2; 3; 4; 5; 5]

# A1 = [1 0 0 0 0 0;
#       1 0 0 0 0 0;
#       0 0 1 0 0 0;
#       0 0 1 0 0 0;
#       0 0 1 0 0 0;
#       0 0 0 0 0 1]
#
# A2 = [1 0 0 0 0 0;
#       0 1 0 0 0 0;
#       0 0 1 0 0 0;
#       0 0 0 1 0 0;
#       0 0 0 0 1 0;
#       0 0 0 0 1 0]
#
# if cutpoint == 7 (second column)
crossover!(R1, R2, support)
