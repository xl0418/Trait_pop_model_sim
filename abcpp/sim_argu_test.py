import sys


def simtest(files,num_threads,result,sstats,iterations,particles):
    orig_stdout = sys.stdout
    f = open('out.txt', 'w')
    sys.stdout = f
    for i in range(100):
        print(i)
    sys.stdout = orig_stdout
    f.close()
    return str(files)+'__' +str(result)+'__'+ str(num_threads)+'__'+str(sstats)+ '__'+str(
        iterations)+'__'+str(particles)