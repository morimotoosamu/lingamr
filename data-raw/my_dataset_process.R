## code to prepare `my_dataset_process` dataset goes here

usethis::use_data(my_dataset_process, overwrite = TRUE)

# import numpy as np
# import pandas as pd
#
# np.random.seed(100)
#
# x3 = np.random.uniform(size=10000)
# x0 = 3.0*x3 + np.random.uniform(size=10000)
# x2 = 6.0*x3 + np.random.uniform(size=10000)
# x1 = 3.0*x0 + 2.0*x2 + np.random.uniform(size=10000)
# x5 = 4.0*x0 + np.random.uniform(size=10000)
# x4 = 8.0*x0 - 1.0*x2 + np.random.uniform(size=10000)
# x6 = 2.0*x1 + np.random.uniform(size=10000)
# x7 = 1.5*x4 + 3.0*x0 + np.random.uniform(size=10000)
# x8 = 0.5*x2 + 2.0*x5 + np.random.uniform(size=10000)
# x9 = 7.0*x3 + 1.0*x6 + np.random.uniform(size=10000)
#
# X = pd.DataFrame(
#   np.array([x0, x1, x2, x3, x4, x5, x6, x7, x8, x9]).T ,
#   columns=['x0', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'x7', 'x8', 'x9']
# )
