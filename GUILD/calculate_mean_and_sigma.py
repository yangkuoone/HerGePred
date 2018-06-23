import math, sys

def main():
    # sc = float(sys.argv[1])
    # list = [ float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]) ]

    # list = map(float, sys.argv[2:])
    sc = 1.2
    list = [1, 2, 3]
    m, s = cal_mean_sigma(list)



def mean(list):
    return sum(list)/len(list)

def sigma(list):
    n = len(list)
    m = mean(list)

    temp = [x*x for x in list]
    msq = sum(temp)/n
    s = msq-m*m
    try:
        s = math.sqrt(s)
    except Exception as e:
        if s < 0.0000000001:
            s = 0
        else:
            raise e
    return s

def cal_mean_sigma(list):
    return mean(list), sigma(list)


if __name__ == "__main__":
    main()

