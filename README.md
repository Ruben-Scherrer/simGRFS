# README 

**simGRFS** is a program to simulate isotropic Gaussian random fields on the sphere with the option to choose from different simulation methods.

For base functionality, no further packages are required. However, for displaying the simulated values on the sphere, the packages *rgl* and *plotrix* are required.

## Basic Usage

There are three main functions:

### 1. SimGRF()

This is the main function to use for simulation.

This function requires as input at least a matrix of points on the sphere. A possible input format is a matrix with columns in format c(x,y,z), where *(x,y,z)* are points on the sphere in cartesian coordinates. A basic example of usage is:

```r
SimGRF(points)
```
This function call leads to a simulation of an isotropic GRF using a default model (see section about GRF models) and the default matrix decomposition method (see section about simulation methods).

Per default, the function returns a matrix that contains 4 rows. Each column has format c(x,y,z,values) where *(x,y,z)* are the points provided, and *values* the simulated GRF values. For example (using 3 random points):

```r
           [,1]       [,2]        [,3]
[1,] -0.3403263 -0.7283710 -0.38687067
[2,] -0.4149311  0.4372996  0.92166460
[3,]  0.8438070 -0.5274891  0.02941866
[4,]  1.5846621 -0.9646873 -0.99784764
```

If only the simulated values of the GRF are needed, the argument *out_format = "values"* can be provided to the function.

### 2. RandomSimGRF()

This function can be used to simulate random points on the sphere. It uses the implemented function *RandomSpherePoints()* to sample randomly on the sphere and then passes these points to *SimGRF()*. The following two chunks do the same.

```r
RandomSimGRF(50)
```

```r
randompoints <- RandomSpherePoints(50)
SimGRF(randompoints)
```

### 3. MultiSimGRF()

This function is used to generate many realizations of the same GRF using optimized methods.

For example, when using the (default) matrix decomposition method, the covariance matrix of the provided points has to be computed only once. This way, multiple realizations can be simulated very fast.

The following code simulates 50 random points 1000 times and returns a matrix, where each row contains the simulated values of all input points.

```r
randompoints <- RandomSpherePoints(50)
MultiSimGRF(1000, randompoints)
```

### GRF models
It is assumed that all GRFs to be simulated have unit variance. An isotropic GRF model is therefore fully specified by its covariance function (a function dependent on $\alpha$, the great circle distance between two points on the sphere).

To specify a model to use for simulation, create an object of the implemented class **GRFModel**.

```r
model1 <- GRFModel()
```

This creates per default an exponential model with covariance $C(\alpha) = \exp(-2\alpha)$. 

Such a model can then be passed as an argument to the simulation function:

```r
SimGRF(points, model1)
```

There is a variety of models implemented to use. For example, the following call creates a GRF model using a $C^2$-Wendland covariance of $C(\alpha)=(1+\rho_1 \rho_2 \alpha)(1-\rho_1 \alpha)_+^{\rho_1}$ with $\rho_1=1, \rho_2=4$.

```r
model2 <- GRFModel(covmodel = "wendc2", param1 = 1, param2 = 4)
```

### GRF model types
There are two **types** of grfmodels:

1. "covar" - utilizes the covariance functions directly

2. "spectral" - utilizes the angular power spectrum (also called Schoenberg coefficients)

The type of a GRF model can be specified by passing "covar" or "spectral" to the argument "modeltype" when creating a model using GRFModel().

```r
model3 <- GRFModel(modeltype = "spectral")
```

The methods "md" and "atbm" need GRF models of type "covar", while the other methods require GRF models of type "spectral". 


### Simulation methods

There are 6 different simulation methods implemented:

1. Matrix Decomposition ("md")

2. Karhunen-Loève Expansion ("kle")

3. Simple Spectral Method ("spectral")

4. Emery Turning Bands Method ("etbm")

5. Different Implementation of Emery Turning Bands Method ("etbm2")

6. Alternative Turning Bands Method ("atbm")

To use different simulation methods, one can pass the requested method as an argument:

```r
SimGRF(points, method = "kle")
```

Note that if the method is chosen manually and a GRF model is provided, one has to take care that the type of the GRF model matches the method.

### Displaying

Instead of returning the simulated values, it is also possible to let the values be displayed on a sphere using colors to indicate the values. This can be done by setting the argument *display=TRUE*, as in the example:

```r
SimGRF(points, method = "kle", display = TRUE)
```

#### Infos
Ruben Scherrer, August 2022




