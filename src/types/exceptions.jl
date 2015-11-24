# This file is part of Kpax3. License is MIT.

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

"""
# Kpax3 Exception

## Description

Exception for wrong data read from a source.

## Fields

* `msg` Optional argument with a descriptive error string
"""
type Kpax3InputError <: Exception
  msg::AbstractString
end
Kpax3InputError() = Kpax3InputError("")
