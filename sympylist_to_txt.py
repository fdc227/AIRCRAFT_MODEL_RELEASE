def sympylist_to_txt(slist, filename):

    file1 = open(filename, '+w')
    m = len(slist)

    try:
        n = len(slist[0])
        if n > 1 and m > 1:
            for i in range(m):
                if i != m - 1:
                    for j in range(n):
                        if j != n - 1:
                            file1.write(str(slist[i][j]))
                            file1.write(',')
                        else:
                            file1.write(str(slist[i][j]))
                            file1.write('\n')
                else:
                    for j in range(n):
                        if j != n - 1:
                            file1.write(str(slist[i][j]))
                            file1.write(',')
                        else:
                            file1.write(str(slist[i][j]))
        elif m > 1:
            for i in range(m):
                if i != m - 1:
                    file1.write(str(slist[i][0]))
                    file1.write('\n')
                else:
                    file1.write(str(slist[i][0]))
    except:
        for j in range(m):
            if j != m - 1:
                file1.write(str(slist[j]))
                file1.write(',')
            else:
                file1.write(str(slist[j]))


    file1.close()       