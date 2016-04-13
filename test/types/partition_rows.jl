# This file is part of Kpax3. License is MIT.

@test_throws KDomainError EwensPitman(1.0, 0.0)
@test_throws KDomainError EwensPitman(2.0, 0.0)

@test_throws KDomainError EwensPitman(0.5, -1.0)
@test_throws KDomainError EwensPitman(0.5, -0.5)

@test_throws KDomainError EwensPitman(0.0, -1.0)
@test_throws KDomainError EwensPitman(0.0, 0.0)

@test_throws KDomainError EwensPitman(0.5, 0)
@test_throws KDomainError EwensPitman(0.0, 0)
@test_throws KDomainError EwensPitman(0.5, 1)
@test_throws KDomainError EwensPitman(0.0, 1)

@test_throws KDomainError EwensPitman(-1.0, 0)
@test_throws KDomainError EwensPitman(-1.0, -1)
@test_throws KDomainError EwensPitman(-1.0, 1.0)
