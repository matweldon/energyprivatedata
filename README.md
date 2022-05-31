Homomorphic regression
================
Mat Weldon
2022-05-31

This project tests whether homomorphic encryption, as implemented in the
`homomorpheR` package, can be used to solve a simplified data sharing
problem.

``` r
library(homomorpheR)
keyPair <- PaillierKeyPair$new(modulusBits = 256)
```

Set up convenience functions to encrypt and decrypt integers

``` r
encrypt = function(x,kp=keyPair){
  if(length(x)>1){
    res = lapply(x,function(y)kp$pubkey$encrypt(y))
  }else{
    res = kp$pubkey$encrypt(x)
  }
  return(res)
}

decrypt = function(x,kp=keyPair){
  if(length(x)>1){
    res = lapply(x,function(y)kp$getPrivateKey()$decrypt(y))
  }else{
      res = kp$getPrivateKey()$decrypt(x)
  }
  return(res)
}
```

Under the Paillier system, multiplying two encrypted quantities is
equivalent to adding the unencrypted quantities.

``` r
a = gmp::as.bigz(12738)
b = gmp::as.bigz(48231)

ea = encrypt(a)
eb = encrypt(b)

identical(decrypt(ea*eb), a+b)
```

    ## [1] TRUE

Additionally, raising an encrypted quantity to some power *n* is
equivalent to multiplying the unencrypted quantity by *n*.

``` r
identical(decrypt(ea^2),a*2)
```

    ## [1] TRUE

This means we can calculate the dot product of an encrypted vector with
an unencrypted vector.

``` r
encrypted_dot_product = function(enc_x,y){
  N = length(enc_x)
  stopifnot(N==length(y))
  
  list_of_products = vector(N,mode="list")
  
  for(i in seq_along(list_of_products)){
    list_of_products[[i]] = enc_x[[i]]^y[i]
  }
  
  dot_prod = Reduce(`*`,list_of_products)
  
  return(dot_prod)
}
```

``` r
x = lapply(1:10, function(x) random.bigz(nBits = 16))
encrypted_x = encrypt(x)
numeric_x = sapply(x,as.integer)
y = 1:10

dot_prod = encrypted_dot_product(encrypted_x,y) |> decrypt()


print(sum(numeric_x*y))
```

    ## [1] 2340793

``` r
print(as.numeric(dot_prod))
```

    ## [1] 2340793

This semi-encrypted dot product would be immediately useful. It would
enable one agency to share encrypted information with another agency,
which could then perform computations and return the result. For
example, agencies could perform regression without sharing data, except
matching id’s.

For example, say ONS has a matrix of sensitive variables about
households, identified by their address, and energy companies have
household energy consumption for those households. ONS could compute and
then elementwise encrypt the matrix:

*ξ*(𝒳) = *ξ*(\[**X**<sup>*T*</sup>**X**\]<sup>−1</sup>**X**<sup>*T*</sup>)

where *ξ*(⋅) represents elementwise encryption. The encryption protects
each element of **X**, but additionally pre-multiplying the design
matrix with the covariance matrix effectively obfuscates the values of
**X**, especially if there are many columns in the matrix.

Then the energy company could compute:

*ξ*(𝒳) ⊙ **y** = *ξ*(\[**X**<sup>*T*</sup>**X**\]<sup>−1</sup>**X**<sup>*T*</sup>**y**) = *ξ*(*β̂*)
where ⊙ represents the homomorphic matrix product derived from the
homomorphic dot product defined above. The energy company would then
send the encrypted estimate vector back to ONS, who would decrypt the
coefficient estimates and combine them with the estimates from other
providers.

This is as far as I’ve been able to get so far because encryption and
decryption can only handle integers or rational numbers.

``` r
a = 1.245

new_a = encrypt(a) |> decrypt()
```
