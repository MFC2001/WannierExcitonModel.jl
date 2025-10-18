# WannierExcitonModel.jl

A Julia package solving Excitonic Bethe-Salpeter equation(BSE) based on Electronic Wannier model.

## Package Features

This package provides a general solution for establishing a BSE model based on electronic Wannier functions, allowing users to customize the model while also implementing several typical scenarios.

## In development

At present, we require users to provide the interaction between wannier bases.
For me, I only know [RESPACK](https://www.sciencedirect.com/science/article/pii/S001046552030391X) give this function, but it's a little bit unstable.
We are trying to calculate the interaction between wannier bases based on well-used open source software, such as [BerkeleyGW](https://berkeleygw.org/), by ourself.
