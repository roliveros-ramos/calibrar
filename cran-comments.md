# version 0.2.0

**round 3**

> Thanks, on CRAN now.
>
> Best
> -k

**round 2**

> Thanks, can you please enclose the URL in
  <http://.....>.

The DESCRIPTION file has been
  updated.

> Best,
Uwe Ligges

Best regards.

**round 1**

> Pls write "This package allows ..."

It's been
  corrected.

> Also, we get
> Undefined global
  functions or variables:
  as.relistable optim pnorm
  read.csv relist rexp > rnorm rpois runif
  weighted.mean write.csv
> Consider adding
>
  importFrom("stats", "optim", "pnorm", "rexp",
  "rnorm", "rpois", "runif", "weighted.mean")
  importFrom("utils", "as.relistable", "read.csv",
  "relist", "write.csv")
to your NAMESPACE
  file.

All those imports have been added to the
  NAMESPACE.

> Pls fix.  Also,
> * checking
  examples ... [81s/81s] OK
> Examples with CPU or
  elapsed time > 5s
                   user system
  elapsed
> calibrar-package 39.012  0.876  39.930
>
  calibrate        36.712  0.968  37.708
> and the
  CRAN Policy asks for only a few seconds: can you pls
  reduce the above run times accordingly?

The running time cannot be reduced, I've added a \dontrun{} to the examples.

> Best
-k

Thanks for the corrections.

Best regards.
Ricardo