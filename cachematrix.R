## Define makeCacheMatrix function

makeCacheMatrix <- function(x = matrix()) {


## x should be a square invertible matrix


## this function returns a list containing functions to:
## 1. set the matrix (set)
## 2. get the matrix (get)
## 3. set the inverse of the matrix (setinver)
## 4. get the inverse of the matrix (getinver)

inver <- NULL ## initialize a vector where inverse will be stored
        set <- function(y) {
                x <<- y
                inver <<- NULL
        }
        get <- function() x
        setinver <- function(inverse) inver <<- mean
        getinver <- function() inver
        list(set = set, get = get, setinver = setinver,getinver = getinver)
}





cacheSolve <- function(x, ...) {
        

## here x is the output of makeCacheMatrix function
## returned should be the inverse of the input matrix to 
## makeCacheMatrix function

inver <- x$getinver()
        
        ## if inverse is already calculated, get it     

        if(!is.null(inver)) {
                message("getting cached data")
                return(inver)
        }

        ## if inverse is not calculated, compute using 'solve' function

        matri.data <- x$get()
        inver <- solve(matri.data, ...)
        x$setinver(inver) ## set the value of inverse
        return(inver) ## Return a matrix that is the inverse of 'x'


}