import argparse
# Runge Kutta Method for solving Equations
# all values returned are in used in 6 decimal places


class SecondOrderRungeKutta:
    """
    This class is a representation of the Second Order Runge Kutta method
    of solving equations in Engineering Mathematics.
    :param eqn: this is the differential equation to be solved by Runge Kutta Method e.g "(x×y’) + y"
                where y’ is the first differential of y (i.e dy/dx)
    :param x: this is the initial state value for `x`
    :param y: this is the initial state value for `y`
    :param y_prime: this is the initial state value for `y'` (i.e first differential of y)
    :param h: this is the step for changing the value of x during the iteration
    """
    iterations = 6  # this determines the number of iterations desired
    decimals = 7

    def __init__(self, eqn, x, y, y_prime, h, upper_limit):
        self.eqn = SecondOrderRungeKutta.__clean_eqn(eqn)
        self.x0 = x
        self.y0 = y
        self.y_prime0 = y_prime
        self.upper_limit = upper_limit
        self.h = h
        self.x = 0
        self.y = 0
        self.y_prime = 0
        self.k1 = 0
        self.k2 = 0
        self.k3 = 0
        self.k4 = 0
        self.p = 0
        self.q = 0
        print(self.eqn)

    def f(self):
        x, y, yprime = (self.x0, self.y0, self.y_prime0) if self.x == 0 else (self.x, self.y, self.y_prime)
        return eval(self.eqn)

    @staticmethod
    def __clean_eqn(eqn):
        operators = {
            '^': '**',
            '{': '(',
            ']': ')',
            '}': ')',
            '[': '(',
            '÷': '/',
            '×': '*',
            "'": 'prime',
            '‘': 'prime',
            '’': 'prime',
            '‛': 'prime'
            }
        for operator in list(operators.keys()):
            eqn = eqn.replace(operator, operators[operator])
        return eqn

    def first_K(self):
        self.k1 = ((self.h**2)/2) * self.f()
        return self.k1

    def second_K(self):
        self.x = self.x0+(self.h/2)
        self.y = round(self.y0+((self.h*self.y_prime0)/2)+(self.k1/4), 6)
        self.y_prime = (self.y_prime0+(self.k1/self.h))
        self.k2 = ((self.h**2)/2) * self.f()
        return self.k2

    def third_K(self):
        self.x = self.x0+(self.h/2)
        self.y = self.y0+((self.h*self.y_prime0)/2)+(self.k1/4)
        self.y_prime = self.y_prime0+(self.k2/self.h)
        self.k3 = ((self.h**2)/2) * self.f()
        return self.k3

    def fourth_K(self):
        self.x = self.x0+self.h
        self.y = self.y0 + (self.h*self.y_prime0) + self.k3
        self.y_prime = (self.y_prime0+((2*self.k3)/self.h))
        self.k4 = ((self.h**2)/2) * self.f()
        return self.k4

    def p_factor(self):
        self.p = (self.k1+self.k2+self.k3)/3
        return self.p

    def q_factor(self):
        self.q = (self.k1+(2*self.k2)+(2*self.k3)+self.k4)/3
        return self.q

    def y1_factor(self):
        y1 = (self.y0+(self.h*self.y_prime0)+self.p)
        return y1

    def y1_prime(self):
        y_prime1 = self.y_prime0 + (self.q/self.h)
        return y_prime1

    def yPPrime(self):
        ypprime = self.f()
        return ypprime

    @classmethod
    def solve(cls, eqn, xNaught, yNaught, yPrime, h, upper_limit=2.0):
        print()
        i = 1
        while xNaught <= upper_limit:
            header = f"ITERATION {i-1}"
            print(f"{header}\n{'=' * len(header)}")
            function = cls(eqn, xNaught, yNaught, yPrime, h, upper_limit)

            print(f"x{i - 1}".ljust(5) + "=\t" + str(round(xNaught, function.decimals)))
            print(f"y{i - 1}".ljust(5) + "=\t" + str(round(yNaught, function.decimals)))
            print(f"y{i - 1}’".ljust(5) + "=\t" + str(round(yPrime, function.decimals)))
            print()
            print(f"y{i-1}’’".ljust(5) + "=\t" + str(round(function.yPPrime(), function.decimals)))
            print()
            print(f"k1".ljust(5) + "=\t" + str(round(function.first_K(), function.decimals)))
            print(f"k2".ljust(5) + "=\t" + str(round(function.second_K(), function.decimals)))
            print(f"k3".ljust(5) + "=\t" + str(round(function.third_K(), function.decimals)))
            print(f"k4".ljust(5) + "=\t" + str(round(function.fourth_K(), function.decimals)))
            print()
            print(f"P".ljust(5) + "=\t" + str(round(function.p_factor(), function.decimals)))
            print(f"Q".ljust(5) + "=\t" + str(round(function.q_factor(), function.decimals)))

            xNaught += function.h
            xNaught = round(xNaught, 2)
            yNaught = function.y1_factor()
            yPrime = function.y1_prime()

            i += 1
            print()
            print()


class FirstOrderRungeKutta:
    pass


def main():
    parser = argparse.ArgumentParser(description="Solve second-order differential equations using the Runge-Kutta method leaving all answers at 6d.p.")
    parser.add_argument("--equation", type=str, help="The differential equation to solve.")
    parser.add_argument("--time_step", "--h", type=float, help="Time step for numerical integration i.e. h")
    parser.add_argument("--xNaught", type=float, help="Naught value for x i.e. value of x at 0")
    parser.add_argument("--yNaught", type=float, help="Naught value for y i.e. value of y at 0")
    parser.add_argument("--yPrime", type=float, help="Value for y prime i.e. value of first differential of y")
    parser.add_argument("--upper_limit", type=float, help="The maximum limit attainable by x after the iterations. Default = 2.0")
    # parser.add_argument("--output", type=str, help="Output file to save results.")

    args = parser.parse_args()
    if args.upper_limit is not None:
        SecondOrderRungeKutta.solve(args.equation, args.xNaught, args.yNaught, args.yPrime, args.time_step, args.upper_limit)
    else:
        SecondOrderRungeKutta.solve(args.equation, args.xNaught, args.yNaught, args.yPrime, args.time_step)

    # Save the result to the output file
    # with open(output_file, 'w') as file:
    #     file.write(result)

if __name__ == "__main__":
    main()