"""part 1 - Q.2"""

def calcDerivative(func):
    """
    :param func: original function
    :return: None
    """
    x = sp.symbols('x')
    f_prime = func.diff(x)
    return f_prime


def simpson(f, startPoint, endPoint, parts):
    """
    :param f: original function
    :param startPoint: start of range
    :param endPoint: end of range
    :param parts: amount of segments
    :return: approximate area of the integral
    """
    if parts % 2 == 1:  # if there is not even numbers of parts
        print("Amount of parts must be even")
        return None
    x = sp.symbols('x')
    func = lambdify(x, f)
    gap = abs(endPoint - startPoint) / parts  # calculate h
    string = "Integral(" + str(startPoint) + ", " + str(endPoint) + ") = 1/3 * " + str(gap) + "[f(" + str(startPoint) + ")"
    appr = func(startPoint)  # placing the start point in the function
    for i in range(1, parts):  # run over the parts
        if i % 2 == 0:  # if is the even place
            string += " + 2 * f(" + str((i * gap) + startPoint) + ")"
            appr += 2 * func((i * gap) + startPoint)
        else:  # if is not the even place
            string += " + 4 * f(" + str((i * gap) + startPoint) + ")"
            appr += 4 * func((i * gap) + startPoint)
        if i % 4 ==0:  # for the printing
            string += "\n"
    string += " * f(" + str(endPoint) + ")]\n"
    print(string)  # print the equation
    appr += func(endPoint)
    appr *= 1 / 3 * gap
    return appr


def rombergMethod(f, a, b, end, epsilon):
    """
    :param f: Original function
    :param a: start of the range
    :param b: end of the range
    :param end: limit of iteration
    :param epsilon: allowed error
    :return: The area in the range
    """
    results = [[0 for i in range(end + 1)] for j in range(end + 1)]  # build matrix
    for k in range(0, end):
        print("R" + str(k+1) + "," + str(1) + " = ", end="")
        res = trapezoidMethod(f, a, b, 2 ** k)  # calculate the values of trapezoid method
        results[k+1][1] = res  # save the value in the matrix
        print(" = " + str(res))  # print the value
    for j in range(2, end + 1):
        for k in range(j, end + 1):
            results[k][j] = results[k][j - 1] + ((1 / ((4 ** (j - 1)) - 1)) * (results[k][j - 1] - results[k - 1][j - 1]))
            print("R" + str(k) + "," + str(j) + " = " + str(results[k][j]))  # print the value
            if abs(results[k][j] - results[k - 1][j]) <= epsilon:  # if the difference is less then epsilon
                return results[k][j]
    return results[k][j]



def trapezoidMethod(f, a, b, n):
    """
    :param f: Original function
    :param a: start of the range
    :param b: end of the range
    :param n: the number of the segments
    :return: The area in the range
    """
    x = sp.symbols('x')
    f = lambdify(x, f)
    h = (b - a) / n
    sum = 0
    save = a
    count = 0
    while a < b:
        sum += 0.5 * ((a + h) - a) * (f(a) + f(a + h))
        count +=1
        if a is not save:
            print(" + ", end="")
        if count is 3:
            print("\n       ", end="")
            count = 0
        print("1/2 * (" + str(b) + " - " + str(a) + ") * (f(" + str(a) + " + f(" + str(a + h) + "))", end="")
        a += h
    return sum


def rangeDivision(polinom, start_point, end_point, epsilon, function):
    """
    :param polinom: Original function
    :param start_point: int value, the start point of the range
    :param end_point: int value, the end point of the range
    :param epsilon: The excepted error
    :param function: The function that needs to be activated
    :return: None
    """
    results = []
    sPoint = start_point
    ePoint = end_point
    flag = False
    temp = start_point + 0.1
    while temp <= end_point:  # while we dont reach to the end point of the range
        print("Range: [" + str(start_point) + ", " + str(temp) + "]")
        res, iter = function(polinom, start_point, temp,
                             epsilon)  # activates the requested function with the original function
        if iter is not None:  # if the return iteration value is not None
            if (res > -epsilon) and (res < epsilon):  # check if the result is very close to 0
                res = 0
            print("The root is " + calcFinalResult(str(res), 10**-4, '13', '18', '41') + "\nNumber of iteration is: " + str(iter))
            results.append(res)
            flag = True
        else:
            print("* There is no Intersection root in range ")
        der = calcDerivative(polinom)  # calculate the derivative
        x = sp.symbols('x')
        res, iter = function(der, start_point, temp,
                             epsilon)  # activates the requested function with the derivative function.
        if iter is not None:  # if the return iteration value is not None
            if (res > -epsilon) and (res < epsilon):  # check if the result is very close to 0
                res = 0
            if (lambdify(x, polinom)(res) > - epsilon) and (
                    lambdify(x, polinom)(res) < epsilon):  # only if the res is root
                print("The root is " + calcFinalResult(str(res), 10**-4, '13', '18', '41') + "\nNumber of iteration is: " + str(iter))
                results.append(res)
                flag = True
            else:
                print("This point is a root of the derivative but is not a root of the function..\n")
        else:
            print("* There is no touching root in range\n")
        start_point = temp  # increase the start point of the range
        temp += 0.1  # increase the end point of the range
    if flag is False:
        print("There is no root in the range " + "[" + str(sPoint) + ", " + str(ePoint) + "]")
    return results


def NewtonRaphson(polinom, start_point, end_point, epsilon):
    """
    :param polinom: Original function
    :param start_point: int value, the start point of the range
    :param end_point: int value, the end point of the range
    :param epsilon: The excepted error
    :return: None
    """
    return rangeDivision(polinom, start_point, end_point, epsilon, calcByNewtonRaphson)


def calcByNewtonRaphson(pol, startPoint, endPoint, epsilon):
    """
    :param pol: Original function
    :param startPoint: int value, the start point of the range
    :param endPoint: int value, the end point of the range
    :param epsilon: The excepted error
    :return: The result and the number of the iteration for getting that result
    """
    Xr = (startPoint + endPoint) / 2  # middle of the range
    iteration = 1
    der = calcDerivative(pol)  # calculate the derivative
    x = sp.symbols('x')
    pol = lambdify(x, pol)
    der = lambdify(x, der)
    if pol(startPoint) * pol(endPoint) > 0:  # check if the is change in the sign in the function
        return None, None
    print("==== Iterations ====")
    while iteration < 100:
        res1 = pol(Xr)
        res2 = der(Xr)
        print("Iteration number: " + str(iteration) + ", Xr = " + str(Xr) + ", f(x) = " + str(res1) + ", f`(x) = " + str(res2))
        iteration += 1
        Xnext = Xr - (res1 / res2)
        if abs(Xnext - Xr) < epsilon:
            print("Iteration number: " + str(iteration) + ", Xr = " + str(Xr) + ", f(x) = " + str(
                res1) + ", f`(x) = " + str(res2))
            print("==== End of iterations ====\n")
            return Xnext, iteration
        Xr = Xnext
    print("The system does not converge... :(")
    return None, None


def secant_method(polinom, start_point, end_point, epsilon):
    """
    :param polinom: Original function
    :param start_point: int value, the start point of the range
    :param end_point: int value, the end point of the range
    :param epsilon: The excepted error
    :return: None
    """
    return rangeDivision(polinom, start_point, end_point, epsilon, calcBySecant)


def calcBySecant(polinom, start_point, end_point, epsilon):
    """
    :param polinom: Original function
    :param start_point: int value, the start point of the range
    :param end_point: int value, the end point of the range
    :param epsilon: The excepted error
    :return: The result and the number of the iteration for getting that result
    """
    Xr = start_point
    Xnext = end_point
    iteration = 1
    x = sp.symbols('x')
    polinom = lambdify(x, polinom)
    if polinom(start_point) * polinom(end_point) > 0:  # check if the is change in the sign in the function
        return None, None
    print("==== Iterations ====")
    while iteration < 100:
        res1 = polinom(Xr)
        res2 = polinom(Xnext)
        print("Iteration number: " + str(iteration) + ", Xr = " + str(Xr) + ", f(x) = " + str(res1) + ", f`(x) = " + str(res2))
        iteration += 1
        temp = Xnext
        Xnext = ((Xr * res2) - (Xnext * res1)) / (res2 - res1)
        Xr = temp
        if abs(Xnext - Xr) < epsilon:
            print("Iteration number: " + str(iteration) + ", Xr = " + str(Xr) + ", f(x) = " + str(
                res1) + ", f`(x) = " + str(res2))
            print("==== End of iterations ====\n")
            return Xnext, iteration
    print("The system does not converge... :(")
    return None, None


def checkDiffer(l, d, epsilon):
    print("**** check the difference between the methods: ****")
    flag = True
    for _ in range(len(l)):
        print("Root " + str(_) + ":\nSecant: " + str(l[_]) + ", Newton Raphson: " + str(d[_]))
        if abs(l[_] - d[_]) > epsilon:
            flag = False
            print("The difference is bigger than the epsilon for some of the roots")
            return
    print("The difference is smaller than the epsilon for all the roots")


def calcFinalResult(result, epsilon, day, hour, minutes):
    """
    :param result: the result
    :param epsilon: the epsilon we decided on for the question
    :param day: the day of the calculation
    :param hour: the hour of the calculation
    :param minutes: the minutes of the calculation
    :return: the result by the requested format
    """
    stringRes = str(result)  # cast the result to string
    i = 0
    while stringRes[i] is not ".":  # run over the string while we get to the point
        i += 1  # count how many digits there is before the point
    i += 1
    count = 1
    while epsilon < 1:  # checking how digit needs after the point
        epsilon *= 10
        count += 1
    differ = len(stringRes)
    while differ < i + count:   # fill difference with zero if the number is to short after the point
        stringRes += "0"
        differ += 1
    stringRes = stringRes[:i + count] + "00000" + day + hour + minutes
    return stringRes


def driver():
    print("******* Part 1 - Q2_a *******")
    x = sp.symbols('x')
    f = (sp.cos((x ** 2) + (5 * x) + 6) / (2 * ((math.exp(1)) ** (-x))))
    print("-----------Secant Method -----------")
    l = secant_method(f, -1.5, 2, 0.0001)
    print("-----------Newton Raphson -----------")
    d = NewtonRaphson(f, -1.5, 2, 0.0001)
    checkDiffer(l, d, 10 ** -7)
    print("\n******* Part 1 - Q2_b *******")

    start_point = 0
    end_point = 1
    part = 18
    epsilon = 10 ** -4
    print("-----Simpson Method -----")
    s = simpson(f, start_point, end_point, part)
    print("Final result:\nIntegral(" + str(start_point) + ", " + str(end_point) + ") = " + str(
        calcFinalResult(s, epsilon, '13', '18', '41')))
    print("\n-----Romberg Method -----")
    r = rombergMethod(f, start_point, end_point, 5, epsilon)
    print("\nFinal result:\nIntegral(" + str(start_point) + ", " + str(end_point) + ") = " + str(
        calcFinalResult(r, epsilon, '13', '18', '41')))

    if abs(s - r) <= epsilon:
        print("\n* The difference between the two methods is smaller than the epsilon")
    else:
        print("\n* The difference between the two methods is bigger than the epsilon - needs another method")



"""part 2 - Q.10"""


def machineEpsilon(func=float):
    machine_epsilon = func(1)
    while func(1)+func(machine_epsilon) != func(1):
        machine_epsilon_last = machine_epsilon
        machine_epsilon = func(machine_epsilon) / func(2)
    return machine_epsilon_last


def driver():
    """
    the main program
    :return: print the results
    """
    print("\n******* Part 2 - Q10_a *******")
    x = sp.symbols('x')
    f = (x * (math.exp(1) ** (-x)) + ln(x ** 2)) * (2 * x ** 3 + 2 * x ** 2 - 3 * x - 5)
    startRange = 0 + machineEpsilon()
    endRange = 1.5
    epsilon = 10 ** (-4)
    print("Newton Raphson method")
    l = NewtonRaphson(f, startRange, endRange, epsilon)
    print("Secant method")
    d = secant_method(f, startRange, endRange, epsilon)
    checkDiffer(l, d, epsilon)
#     ------------ integral ---------------------
    print("\n******* Part 2 - Q10_b *******")
    startRange = 0.5
    endRange = 1
    print("simpson method")
    s = simpson(f, startRange, endRange, 6)
    print("Final result: ", end="")
    print(calcFinalResult(s, epsilon, '13', '18', '52'))
    print("romberg Method")
    r = rombergMethod(f, startRange, endRange, 5, epsilon)
    print("Final result: ", end="")
    print(calcFinalResult(r, epsilon, '13', '18', '33'))

    if abs(s - r) <= epsilon:
        print("\n* The difference between the two methods is smaller than the epsilon")
    else:
        print("\n* The difference between the two methods is bigger than the epsilon - needs another method")

"""part 2 - Q.13"""


def driver():
    """
    the main program
    :return: print the results
    """
    print("\n******* Part 2 - Q13_a *******")
    x = sp.symbols('x')
    f = (2 * x * (math.exp(1) ** (-x)) + ln(2 * x ** 2)) * (2 * x ** 2 - 3 * x - 5)
    startRange = 0 + machineEpsilon()
    endRange = 3
    epsilon = 10 ** (-4)
    print("Newton Raphson method")
    l = NewtonRaphson(f, startRange, endRange, epsilon)
    print("Secant method")
    d = secant_method(f, startRange, endRange, epsilon)
    checkDiffer(l, d, epsilon)
#     ------------ integral ---------------------
    print("\n******* Part 2 - Q13_b *******")
    startRange = 0.5
    endRange = 1
    print("\nsimpson method")
    s = simpson(f, startRange, endRange, 6)
    print("Final result: " + calcFinalResult(s, epsilon, '13', '18', '33'))
    print("\nromberg Method")
    r = rombergMethod(f, startRange, endRange, 5, epsilon)
    print("Final result: " + calcFinalResult(r, epsilon, '13', '18', '33'))

    if abs(s - r) <= epsilon:
        print("\n* The difference between the two methods is smaller than the epsilon")
    else:
        print("\n* The difference between the two methods is bigger than the epsilon - needs another method")



"""part 3 - Q.21"""

def GausssElimination(A, b):
    """
    :param A: the coefficients matrix
    :param b: the solution vector
    :return: none
    """
    print("** Gauss Elimination **")
    tempMatrix = Invers(A)  # A^(-1)
    condA = infinityNorma(tempMatrix)  # calculate the infinity norma of A^(-1)
    printMatrix(tempMatrix, "inverse A")  # print the inverse A
    tempMatrix = multMat(tempMatrix, b)
    printMatrix(tempMatrix, "x")  # print x vector
    condA *= infinityNorma(A)  # calculate the infinity norma of A and find condA
    print("-----------------------------")
    print("condA = " + str(condA))
    print("-----------------------------")
    return tempMatrix


def infinityNorma(mat):
    """
    :param mat: a list - matrix
    :return: the infinity norma
    """
    maximum = 0
    for i in range(len(mat)):  # rum over the lines
        newSum = lineSum(mat[i])  # send to function that calculate the line sum
        maximum = max(maximum, newSum)  # choose the maximum
    return maximum


def lineSum(line):
    """
    :param line: A list od numbers - line for the matrix
    :return: the sum of all the numbers in abs  in the list
    """
    lineSum = 0
    for index in range(len(line)):  # run over all the line`s members
        lineSum += abs(line[index])
    return lineSum


def Invers(A):
    """
    :param A: a list - matrix
    :return: the invers of the matrix
    """
    if det(A) is 0:
        print("A is singular")
        return
    inverseMat, counter = toUpperWithPivoting(A)  # return the multiplication of the elementary matrices that create U
    printMatrix(inverseMat, "E" + str(counter))
    counter += 1
    tempMatrix = multMat(inverseMat, A)  # upper triangle matrix
    mat, c = FromUpperToInvers(tempMatrix, inverseMat, counter)  # return inverse matrix
    return mat


def det(A):
    """
    :param A: a square matrix
    :return: the determinant of A
    """
    if len(A[0]) is not len(A):  # check if A is a square matrix
        return None
    if len(A) == 2:  # if the matrix size is 2*2
        return (A[0][0] * A[1][1]) - (A[0][1] * A[1][0])
    sum = 0  # will save the determinant
    for j in range(len(A)):
        sum += (-1) ** (0 + j) * A[0][j] * det(minor(A, 0, j))
    return sum


def minor(A, i, j):
    """
    :param A: a square matrix
    :param i: a row index for the minor
    :param j: a column index for the minor
    :return: the minor that is created from deleting row i and column j from A
    """
    matSize = len(A)
    mat = newMat(matSize - 1, matSize - 1)  # the mat for the minor
    iMinor = 0
    jMinor = 0
    for iIndex in range(matSize):  # go over the rows of the matrix
        for jIndex in range(matSize):  # go over the columns of the matrix
            if iIndex == i:  # don't copy the row we want to delete
                iMinor -= 1
                break
            elif jIndex != j:  # don't copy the row we want to delete
                mat[iMinor][jMinor] = A[iIndex][jIndex]  # copy the rest of the elements in the matrix
                jMinor += 1
        jMinor = 0
        iMinor += 1
    return mat


def toUpperWithPivoting(A):
    """
    :param A: the matrix we want to turn into a upper triangle matrics
    :return: the multiplication of the elementary matrics that create U
    """
    InL = createIdentityMatrix(len(A))  # creating a identity matrix
    counter = 1
    for i in range(len(A) - 1):  # run over the columns
        maxPivot = abs(A[i][i])
        lineIndex = i
        for j in range(i + 1, len(A)):  # run over the elements undet the pivot
            if abs(A[j][i]) > maxPivot:  # if the element is greater then the pivot
                maxPivot = abs(A[j][i])
                lineIndex = j
        elementary = linesExchange(A, i, lineIndex)
        InL = multMat(elementary, InL)
        printMatrix(elementary, "E" + str(counter))
        counter += 1
        A = multMat(elementary, A)
        if A[i][i] != 0:  # check if B is regular
            for j in range(i + 1, len(A)):  # run over the lines
                identity = createIdentityMatrix(len(A))  # creating identity matrix
                identity[j][i] = -(A[j][i] / A[i][i])  # elementary matrix
                printMatrix(identity, "E" + str(counter))
                counter += 1
                InL = multMat(identity, InL)  # L^(-1)
                A = multMat(identity, A)
    return InL, counter


def driver():
    """
    main function for iterative isolation of variables method
    :return: prints results
    """
    A = [[10, 8, 1],
         [4, 10, -5],
         [5, 1, 10]]

    b = [[-7],
         [2],
         [1.5]]

    epsilon = 10 ** (-4)

    flagDom = dominantDiagonal(A)  # check if there is a dominant diagonal
    name1 = "LUdecomposition"
    name2 = "iterativGuassSeidel"
    name3 = "GausssElimination"
    # check
    if flagDom is False:
        copyA = copyMat(A)
        copyB = copyMat(b)
        copyA, copyB = createDominantDiagonal(copyA, copyB)  # change the matrix to be with dominant diagonal
        if (copyA is not None) and (copyB is not None):  # check if the return matrices are not none
            A = copyA
            b = copyB
            flagDom = True

        # end check
    x1 = LUdecomposition(A, b)
    x2 = iterativGuassSeidel(A, b, epsilon, flagDom)
    x3 = GausssElimination(A, b)
    checkDifferPart3(x1, x2, name1, name2, epsilon)
    checkDifferPart3(x1, x3, name1, name3, epsilon)
    calcFinalResultVector(x1, name1, epsilon, "14", "00", "40")
    calcFinalResultVector(x2, name2, epsilon, "14", "00", "41")
    calcFinalResultVector(x3, name3, epsilon, "14", "00", "42")

"""part 3 - Q.29"""

def LUdecomposition(mat, b):
    """
    :param mat: the coefficients matrix
    :param b: the solution vector
    :return: none
    """
    inversL, L, counter = toUpperTriangularMat(mat)  # calculate the L matrix and the inverse L
    print("-----------------------------")
    string = "L = "
    for i in range(1, counter):
        if i < counter - 1:
            string += "invJ" + str(i) + " * "
        else:
            string += "invJ" + str(i)
    print(string)
    printMatrix(L, "L")  # print the L matrix
    print("-----------------------------")
    string = "inverse L = "
    for i in range(counter - 1, 0, -1):
        if i > 1:
            string += "J" + str(i) + " * "
        else:
            string += "J" + str(i)
    print(string)
    printMatrix(inversL, "inverse L")
    print("-----------------------------")
    print("U = invL * A")
    U = multMat(inversL, mat)  # calculate the U matrix
    inversU, counter = FromUpperToInvers(U, createIdentityMatrix(len(U)))  # calculate thr inverse of U
    print("-----------------------------")
    string = "inverse U = "
    for i in range(counter - 1, 0, -1):
        string += "E" + str(i) + " * "
    string += "U"
    print(string)
    printMatrix(U, "U")  # print the U matrix
    print("-----------------------------")
    print("x = invU * invL * b")
    x = multMat(inversL, b)  # finding the result vector
    x = multMat(inversU, x)  # finding the result vector
    printMatrix(x, "x")  # print the X matrix
    return x


def FromUpperToInvers(A, inverseMat, counter=1):
    """
    :param A: upper matrix
    :param inverseMat: the matrix that will become the inverse
    :return: Inverse matrix
    """
    elemntarMat = createIdentityMatrix(len(A))  # identity matrix
    for i in range(len(A) - 1, -1, -1):  # run over the columns
        for j in range(i):  # run over the lines above the pivot
            elemntarMat[j][i] = -(A[j][i] / A[i][i])
            A = multMat(elemntarMat, A)
            inverseMat = multMat(elemntarMat, inverseMat)
            printMatrix(elemntarMat, "E" + str(counter))
            counter += 1
            elemntarMat[j][i] = 0
        if A[i][i] != 1:  # convert the pivots to one
            elemntarMat[i][i] = 1 / A[i][i]
            A = multMat(elemntarMat, A)
            inverseMat = multMat(elemntarMat, inverseMat)
            printMatrix(elemntarMat, "E" + str(counter))
            counter += 1
            elemntarMat[i][i] = 1
    return inverseMat, counter


def toUpperTriangularMat(A):
    """
    :param A: the matrix we want to turn into a upper triangle matrics
    :return: the multiplication of the elementary matrics that create U
    """
    L = createIdentityMatrix(len(A))  # create indetity matrix
    InL = createIdentityMatrix(len(A))  # create indetity matrix
    counter = 1
    for i in range(len(A)):  # run over the lines
        if A[i][i] == 0:  # if the pivot is 0
            for j in range(i + 1, len(A)):  # run over the columns
                if A[j][i] != 0:  # if the element under the pivot is not 0
                    elementary = linesExchange(A, i, j)
                    L = multMat(L, elementary)  # make lines exchange and multiply
                    printMatrix(elementary, "J" + str(counter))
                    InL = multMat((linesExchange(A, i, j)), InL)  # make lines exchange and multiply
                    printMatrix(identity, "inverse J" + str(counter))
                    A = multMat((linesExchange(A, i, j)), A)
                    counter += 1
                    break
        if A[i][i] != 0:  # check if B is regular
            for j in range(i + 1, len(A)):  # run over the columns
                identity = createIdentityMatrix(len(A))
                identity[j][i] = -(A[j][i] / A[i][i])  # elementary matrix
                printMatrix(identity, "J" + str(counter))
                InL = multMat(identity, InL)  # L^(-1)
                A = multMat(identity, A)
                identity[j][i] *= -1  # changing the element in order to find L
                printMatrix(identity, "inverse J" + str(counter))
                L = multMat(L, identity)
                counter += 1
    return InL, L, counter


def linesExchange(A, line1, line2):
    """
    :param A: A matrix
    :param line1: A line
    :param line2: The line we want to exchange with
    :return: elementry matrix
    """
    idendityMax = createIdentityMatrix(len(A))  # create identity matrix

    # exchange the members in line1
    temp = idendityMax[line1][line1]
    idendityMax[line1][line1] = idendityMax[line2][line1]
    idendityMax[line2][line1] = temp

    # exchange the members in line2
    temp = idendityMax[line2][line2]
    idendityMax[line2][line2] = idendityMax[line1][line2]
    idendityMax[line1][line2] = temp
    return idendityMax


def createIdentityMatrix(size):
    """
    :param size: the size of the square matrix
    :return:
    """
    identityMat = newMat(size, size)  # create a zero matrix in the required size
    for index in range(size):  # go over the main diagonal
        identityMat[index][index] = 1  # change the elements in the main diagonal to 1
    return identityMat


def multMat(A, B):
    """
    :param A: a matrix in sise n*m
    :param B: a mtrix in size m*k
    :return: A*B  (in size n*k)
    """
    if len(A[1]) == len(B):  # check if A and B have the same number of rows and columns
        C = newMat(len(A), len(B[0]))  # the matrix the function returns
        for i in range(len(C)):
            for j in range(len(C[1])):
                for k in range(len(B)):
                    C[i][j] += A[i][k] * B[k][j]
        return C
    else:
        return None  # the multiplication  is impossible


# ############### 1st method ############### #

##############################################

# ############### 2nd method ############### #
def iterativGuassSeidel(A, b, epsilon, flagD):
    """
    :param A: a matrix
    :param b: the result vector
    :param epsilon: the mistake
    :param flagD: tell us if the system have dominant diagonal
    :return: vector x if found
    """
    print("** Iterative Guass Seidel **")
    flagC = False  # flagC = false if the linear equations does not converge
    x = newMat(len(b), 1)  # create the result vector
    print("The results are:\nThe first guess is: ")
    printMatrix(x)
    for _ in range(99):  # max number of iterations is 99
        oldX1 = x[0][0]  # copy the old x value of the current iteration
        for i in range(len(x)):  # go over the all variables
            if A[i][i] == 0:  # preventing division by zero
                return None
            temp = b[i][0] / A[i][i]
            for j in range(len(x)):  # calculate the value of the variable in the new iteration
                if i != j:
                    temp -= (A[i][j] * x[j][0]) / A[i][i]
            x[i][0] = temp  # update the calculated value
        print("The result of the " + str(_ + 1) + " iteration is: ")
        printMatrix(x)
        if abs(oldX1 - x[0][0]) < epsilon:  # check stop condition
            flagC = True
            break
    if flagC is True:
        print("***********")  # print final results
        if flagD is False:
            print("Although there is no dominant diagonal, the results are: \n")
        print("The final result is: x = ")
        printMatrix(x)
        print("The number of the iteration is : " + str(_ + 1))
        return x
    else:
        print("The system of linear equations does not converge :(")


# ############### 2nd method ############### #


def printMatrix(a, name=None):
    """
    :param a: a matrix to print
    :return: prints in matrix format
    """
    print("-----------------------------")
    if name is not None:
        print(name + " = ")
    for i in range(len(a)):
        if i is len(a) - 1:
            print(" " + str(a[i]) + "]")
        elif i is 0:
            print("[" + str(a[i]))
        else:
            print(" " + str(a[i]))
    print("")

def newMat(numRow, numCol):
    """
    :param numRow: the number of rows in the mat
    :param numCol: the number of columns in the mat
    :return: a zero matrix in the required size
    """
    mat = []  # the zero matrix the function returns
    for i in range(numRow):
        mat.append([])  # create a new row
        for j in range(numCol):
            mat[i].append(0)  # fill the row with
    return mat


def dominantDiagonal(A):
    """
    :param A: a list - matrix
    :return: true if the matrix have dominant diagonal else return false
    """
    for i in range(len(A)):
        lineSum = 0  # sum of every line except to the element in the diagonal
        for j in range(len(A)):
            if i != j:
                lineSum += A[i][j]
        if A[i][i] <= lineSum:
            # print("There is no dominant diagonal ...")
            return False
    print("There is a dominant diagonal :)")
    return True


# dominant diagonal part

def copyMat(A):
    """
    :param A: a matrix
    :return: a copy of the matrix
    """
    B = newMat(len(A), len(A[0]))  # create a zero matrix of the same size as A
    for i in range(len(A)):  # copy A
        for j in range(len(A[0])):
            B[i][j] = A[i][j]
    return B


def createDominantDiagonal(A, b=None):
    """
    :param A: a coefficients matrix
    :param b: the column vector of constant terms.
    :return: matrix A with dominant diagonal
    """
    max = 0  # the max value in the current row or column in the matrix
    maxIndex = 0  # the index of the max value
    for i in range((len(A))):  # calc the max value for each member on the main diagonal
        sum = 0  # the sum of the members in the current row in A
        for j in range(len(A)):  # go over the current row
            sum += abs(A[i][j])  # add the value of each member in the row
            if abs(A[i][j]) > max:  # search for the max value in the current row
                max = abs(A[i][j])
                maxIndex = j
        if (sum - max) <= max:  # if the max value in the row meets the condition of a dominant diagonal
            A = manualSwapCol(A, maxIndex,
                              i)  # swap between the columns of the current value on the main diagonal and the max value in that row
        else:  # look for the max value in the current column
            max = 0
            maxIndex = 0
            for j in range(len(A)):  # go over the current column
                # sum += abs(A[j][i])
                if abs(A[j][i]) > max:  # search for the max value in the current column
                    max = abs(A[j][i])
                    maxIndex = j
            if rowSum(A[j]) - max <= max:  # if the max value in the row meets the condition of a dominant diagonal
                A, b = manualSwapRow(A, b, i,
                                     maxIndex)  # swap between the rows of the current value on the main diagonal and the max value in that column
            else:
                print("ERROR - no dominant diagonal")  # A can't be changed into dominant diagonal matrix
                return None, None
    return A, b


def manualSwapRow(a, b, r1, r2):
    """
    manaul rows exchange (without e)
    :param a: The coefficient matrix
    :param b:  The column vector of constant terms
    :param r1: the first row to swap
    :param r2: the second row to swap
    :return: the matrix after the swap, The column vector of constant terms after swap
    """

    if r2 < len(a) and r1 < len(a):
        temp = a[r1]
        a[r1] = a[r2]
        a[r2] = temp
        if b is not None:  # if the result vector is not none swap him too
            temp = b[r1]
            b[r1] = b[r2]
            b[r2] = temp
    return a, b


def manualSwapCol(a, c1, c2):
    """
    :param a: The coefficient matrix
    :param c1: the first column to swap
    :param c2: the second column to swap
    :return: the matrix after the swap
    """
    if c2 < len(a) and c1 < len(a):
        for i in range(len(a)):
            temp = a[i][c1]
            a[i][c1] = a[i][c2]
            a[i][c2] = temp
    return a


def rowSum(line):
    """
    :param line: A list od numbers - line for the matrix
    :return: the sum of all the numbers in abs  in the list
    """
    lineSum = 0
    for index in range(len(line)):  # run over all the line`s members
        lineSum += abs(line[index])
    return lineSum


# end dominant part


def checkDifferPart3(l, d, name1, name2, epsilon):
    print("check the difference between the methods:")
    flag = True
    for i in range(len(l)):
        print("x" + str(i + 1) + " - " + name1 + ": " + str(l[i][0]) + "\nx" +
              str(i + 1) + " - " + name2 + ": " + str(d[i][0]))
        if abs(l[i][0] - d[i][0]) > epsilon:
            flag = False
            print("The difference is bigger than epsilon for some of the components\n")
            return
    print("The difference is smaller than epsilon for all the components\n")



def calcFinalResultVector(vector, name, epsilon, day, hour, minutes):
    """
    :param vector: vector x
    :param epsilon: allowed error
    :param day: dd
    :param hour: hh
    :param minutes: mm
    :return: prints the formatted result vector x
    """
    for i in range(len(vector)):
        vector[i][0] = calcFinalResult(vector[i][0], epsilon, day, hour, minutes)
    printMatrix(vector, name + " - Final Result in Format")


# ############## main ###############
def driver():
    """
    main function for iterative isolation of variables method
    :return: prints results
    """
    A = [[1, 0, -1],
         [-0.5, 1, -0.25],
         [1, -0.5, 1]]

    b = [[0.2],
         [-1.425],
         [2]]

    epsilon = 10 ** (-4)

    flagDom = dominantDiagonal(A)  # check if there is a dominant diagonal
    name1 = "LUdecomposition"
    name2 = "iterativGuassSeidel"
    # check
    if flagDom is False:
        copyA = copyMat(A)
        copyB = copyMat(b)
        copyA, copyB = createDominantDiagonal(copyA, copyB)  # change the matrix to be with dominant diagonal
        if (copyA is not None) and (copyB is not None):  # check if the return matrices are not none
            A = copyA
            b = copyB
            flagDom = True

        # end check
    x1 = LUdecomposition(A, b)
    x2 = iterativGuassSeidel(A, b, epsilon, flagDom)
    checkDifferPart3(x1, x2, name1, name2, epsilon)
    calcFinalResultVector(x1, name1, epsilon, "14", "00", "55")
    calcFinalResultVector(x2, name2, epsilon, "14", "00", "56")




"""part 4 - Q.35"""

def neville(pointsList, m, n, X):
    """
    :param pointsList: a list of points
    :param m: the start point index
    :param n: the end point index
    :param X: an x value for it we want to find the y value in the function that is created from the given points
    :return: a proximity of f(x)
    """
    if m is n:  # the start index and end index are the same (the same point)
        print("P" + str(m) + " = " + str(pointsList[m][1]))
        return pointsList[m][1]   # return the y value of the point
    else:
        Xm = pointsList[m][0]  # the start point x value
        Xn = pointsList[n][0]  # the end point x value
        res = (X - Xm) * neville(pointsList, m+1, n, X) - (X - Xn) * neville(pointsList, m, n-1, X)
        res /= (Xn - Xm)
        print("P" + str(m) +"," + str(n) + " = " + str(res))
        return res   # the result (f(X))


def polynomial(pointsList, X):
    """
    :param pointsList: the list of the points
    :param X: the point that we want to find her approximate value
    :return: the approximate value of X
    """
    mat = newMat(len(pointsList), len(pointsList))
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            mat[i][j] = pow(pointsList[i][0], j)  # The coefficient matrix
    print("A = ")
    printMatrix(mat)
    b = newMat(len(pointsList), 1)
    for i in range(len(b)):
        b[i][0] = pointsList[i][1]  # The column vector of constant terms
    print("b = ")
    printMatrix(b)
    matRes = LUdecomposition(mat, b)  # returns the solution matrix
    print("p" + str(len(pointsList)) + "(x) = ", end="")
    string = ""
    for i in range(len(matRes)):
        string += str(matRes[i][0])
        if i is not 0:
            string += " * x^" + str(i)
        if i is not len(matRes) - 1:
            string += " + "
    print(string + "\n")
    # calc mat
    res = 0
    for i in range(len(matRes)):
        res += matRes[i][0] * pow(X, i)  # calc the y value for the requested x
    return res

def Driver():
    point_list = [[1.2, 3.5095], [1.3, 3.6984], [1.4, 3.9043], [1.5, 4.1293], [1.6, 4.3756]]
    X = 1.37
    epsilon = 10 ** -4
    print("==== Neville Method ====")
    n = neville(point_list, 0, len(point_list) - 1, X)
    print("Final result:\nf(" + str(X) + ") = " + str(calcFinalResult(n, epsilon, '13', '18', '59')))
    print("\n==== Polynomial Method ====\n")
    p = polynomial(point_list, X)
    print("Final result:\nf(" + str(X) + ") = " + str(calcFinalResult(p,epsilon ,'13', '18', '59')))

    if abs(n - p) < epsilon:
        print("\n* The difference between the two methods is smaller than the epsilon")
    else:
        print("\n* The difference between the two methods is bigger than the epsilon - needs another method")
