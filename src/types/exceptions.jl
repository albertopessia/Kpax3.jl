# This file is part of Kpax3. License is MIT.

"""
# Kpax3 Exception

## Description

Provides a message explaining the reason of the DomainError exception.

## Fields

* `msg` Optional argument with a descriptive error string
"""
type KDomainError <: Exception
  msg::AbstractString
end
KDomainError() = KDomainError("")

"""
# Kpax3 Exception

## Description

Exception for a wrong formatted FASTA file.

## Fields

* `msg` Optional argument with a descriptive error string
"""
type KFASTAError <: Exception
  msg::AbstractString
end
KFASTAError() = KFASTAError("")

"""
# Kpax3 Exception

## Description

Exception for wrong data read from a source.

## Fields

* `msg` Optional argument with a descriptive error string
"""
type KInputError <: Exception
  msg::AbstractString
end
KInputError() = KInputError("")
