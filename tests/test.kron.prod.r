

# test that for kron.prod
# =========

n1 <- 2
n2 <- 3
n3 <- 3
m1 <- matrix(rnorm(n1^2),nrow=n1,ncol=n1)
m2 <- matrix(rnorm(n2^2),nrow=n2,ncol=n2)
m3 <- matrix(rnorm(n3^2),nrow=n3,ncol=n3)
y <- runif(n=n1*n2*n3)

# build tensor product matrix
make.kron  <- system.time(kron      <- kronecker(m3,kronecker(m2,m1)))

# notice: in kron the indexing corresponds to

expand.grid( fastest=1:nrow(m1), middle=1:nrow(m2), slowest=1:nrow(m3) )


# do multiplication
calc.kron  <- system.time(R.kron <- kron %*% y)
# total time
total.kron <- make.kron + calc.kron
print(total.kron)
# size of tensor product matrix
print(object.size(kron),units="auto")

# run and time kron.prod
print(system.time(my.kron <- kron.prod(y,list(m1,m2,m3))))

# notice again INDEXING!: in kron.prod need to feed the matrices in reverse order as in kron!

# check result
sum(abs(R.kron-my.kron))

# more dimensions: 6
n1 <- 2
n2 <- 5
n3 <- 4
n4 <- 3
n5 <- 3
n6 <- 3
m1 <- matrix(rnorm(n1^2),nrow=n1,ncol=n1)
m2 <- matrix(rnorm(n2^2),nrow=n2,ncol=n2)
m3 <- matrix(rnorm(n3^2),nrow=n3,ncol=n3)
m4 <- matrix(rnorm(n4^2),nrow=n4,ncol=n4)
m5 <- matrix(rnorm(n5^2),nrow=n5,ncol=n5)
m6 <- matrix(rnorm(n6^2),nrow=n6,ncol=n6)
y <- runif(n=n1*n2*n3*n4*n5*n6)

# build tensor product matrix
make.kron  <- system.time(kron      <- kronecker(m6,kronecker(m5,kronecker(m4,kronecker(m3,kronecker(m2,m1))))))
# do multiplication
calc.kron  <- system.time(R.kron <- kron %*% y)
# total time
total.kron <- make.kron + calc.kron
print(total.kron)
# size of tensor product matrix
print(object.size(kron),units="auto")

# run and time kron.prod
print(system.time(my.kron <- kron.prod(y,list(m1,m2,m3,m4,m5,m6))))

# check result
sum(abs(R.kron-my.kron))

# bigger?

n1 <- 30
n2 <- 20
n3 <- 20
m1 <- matrix(rnorm(n1^2),nrow=n1,ncol=n1)
m2 <- matrix(rnorm(n2^2),nrow=n2,ncol=n2)
m3 <- matrix(rnorm(n3^2),nrow=n3,ncol=n3)
y <- runif(n=n1*n2*n3)

# build tensor product matrix
# my computer dies when attempting to build this tensor product

# run and time kron.prod
print(system.time(my.kron <- kron.prod(y,list(m1,m2,m3))))


