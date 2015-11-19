"""
# Kpax3 Exception

## Description

Provides a message explaining the reason of the DomainError exception.

## Fields

* `msg` Optional argument with a descriptive error string
"""
type Kpax3DomainError <: Exception
  msg::AbstractString
end
Kpax3DomainError() = Kpax3DomainError("")

"""
# Kpax3 Exception

## Description

Exception for a wrong formatted FASTA file.

## Fields

* `msg` Optional argument with a descriptive error string
"""
type Kpax3FASTAError <: Exception
  msg::AbstractString
end
Kpax3FASTAError() = Kpax3FASTAError("")
