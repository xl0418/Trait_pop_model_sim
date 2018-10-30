import numpy as np
import dvtraitsim_cpp as dvcpp


hist = np.zeros(20, np.int)
for i in range(0,10000):
    hist[dvcpp.split_binomial50(10)] += 1;
for x in hist:
    print(x)

