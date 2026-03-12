# CLASS 6 R FUNCTIONS!
Karolina A19106745

## Background

All functions in R have at least 3 things:

- A **name** that we use to call the function
- One or more input **arguments**
- The **body** the lines of R code that do the work

## Our first function

Let’s write a silly wee function to `add` some numbers(the input
arguments)

``` r
add <- function(x, y){
x + y
}
```

Now we can use this function

``` r
add (100, 1)
```

    [1] 101

``` r
add(x=c(100, 1, 100), y=1)
```

    [1] 101   2 101

> Q. What if I give a multiple element vector to `x1 and y`

``` r
add(x=c(100, 1), y=c(100, 1))
```

    [1] 200   2

> Q. What if I give three inputs to the function?

``` r
#add(x=c(100, 1), y=1, z=1)
```

> Q. What if I give only one input to the add function?

``` r
addnew<- function(x, y=1){
x + y
}
```

``` r
addnew(x=100)
```

    [1] 101

``` r
addnew(c(100,1), 100)
```

    [1] 200 101

If we write our function with input arguments having default value then
the user does not have to use a user specified value.

## A second function

Let’s try something more interesting: Make a sequence generating tool..

The `sample()` function can be a useful starting point here:

``` r
sample(1:10, size=4)
```

    [1] 9 8 5 2

> Q. Generate 9 random numbers taken from the input vector x=1:10?

``` r
sample(1:10, size=9)
```

    [1]  6  5  2  3 10  4  9  1  8

> Q. Generate 12 random numbers taken from the input vector x=1:10?

``` r
#sample(1:10, size=12)
```

``` r
sample(1:10, size=12, replace=TRUE)
```

     [1] 8 9 6 7 9 6 3 2 2 7 9 6

> Q. Write code for the `sample()` function that generates nucleotide
> sequences of length 6?

``` r
sample(x=c("A", "G", "C", "T"), size=6, replace=TRUE)
```

    [1] "C" "C" "G" "C" "A" "C"

> Q. Write a first function `generate_dna()` that returns a user
> specified length DNA sequence:

``` r
generate_dna <- function(length) {
  
  sample(x=c("A", "G", "C", "T"), size=length, replace=TRUE)
  
}
```

``` r
generate_dna(12)
```

     [1] "A" "G" "A" "G" "C" "A" "T" "T" "C" "T" "A" "C"

> **Key-Points** Every function in R looks fundamentally the same in
> terms of its structure. Basically 3 things: name, input, and body


    name <- function(input){
      
      body
    }

> Functions can have multiple inputs. These can be **required** or
> **optional** arguments. With optional arguments having a set default
> value.

> Q. Modify and improve our `generate_dna()` function to return it’s
> generated sequence in a more standard formt like “AGTAGTA” rather than
> the vector “A”, “C”, “G”, “A”

``` r
generate_dna <- function(length=6, fasta=TRUE) {
  
ans <- sample(x=c("A", "G", "C", "T"), size=length, replace=TRUE)


if(fasta) {
  cat("Single-element vector output")
  ans <- paste(ans, collapse = "")

} else {
  
  cat("Multi-element vector output")
}

  return(ans)  
}

generate_dna()
```

    Single-element vector output

    [1] "ACCCCC"

The `paste()` function - it’s job is to join up or stick together (a.k.a
paste) input strings together

``` r
paste("alice", "loves R", sep=" ")
```

    [1] "alice loves R"

Flow control means where the R brain goes in your code

``` r
good_mood <- TRUE

if(good_mood) {
  
  cat("Great!")
  
  
} else {
  
  cat("Bummer!")
  
}
```

    Great!

``` r
good_mood <- FALSE

if(good_mood) {
  
  cat("Great!")
  
  
} else {
  
  cat("Bummer!")
  
}
```

    Bummer!

## A Protein generating function

> Q. Write a function, called ’generate_protein()\`, that generates a
> user specified length protein sequence

There are 20 natural amino-acids

``` r
aa <- c("A", "R", "N", "D", "B", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
```

``` r
generate_protein <- function(length) {
  # The amino-acids to sample from 
  aa <- c("A", "R", "N", "D", "B", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

  # Draw n-length amino acids to make our sequence
 ans <-  sample(aa, size=length, replace=T)
  ans <- paste(ans, collapse = "")
}
```

``` r
myseq<- generate_protein(42)
myseq
```

    [1] "WQTWDMFRKEIIGHPNLQFPMQYPYPENRTAVITLCHFSBLA"

> Q. Use that function to generate random protein sequences between
> length 6 and 12

``` r
generate_protein(6)
generate_protein(7)
generate_protein(8)
generate_protein(9)
generate_protein(10)
generate_protein(11)
generate_protein(12)
```

``` r
for(i in 6:12) {
  
  # FASTA ID line ">id"
  
  cat(">", i, sep="", "\n")
  
  # Our protein sequence line
  cat(generate_protein(i), "\n")
  
}
```

    >6
    HLGNBR 
    >7
    SDLFDWY 
    >8
    IEADKCWV 
    >9
    TRCCMIMKM 
    >10
    YHHVFDACKB 
    >11
    DFILLNVKNIV 
    >12
    EAPMGLPHTPQR 

> Q. Are any of your sequences unique i.e not found anywhere in nature?

Yes! Amino acids 6 through 12.
