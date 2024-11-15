using System;

class Program
{
    // задана система рівнянь
    static double[] SystemOfEquations(double x, double y)
    {
        return new double[] {
            Math.Tan(x * y + 0.1) - x * x,
            x * x + 2 * y * y - 1
        };
    }

    // матриця Якобі
    static double[,] JacobiMatrix(double x, double y)
    {
        double sec2_xy_01 = 1 / (Math.Cos(x * y + 0.1) * Math.Cos(x * y + 0.1));

        return new double[,] {
            { sec2_xy_01 * y - 2 * x, sec2_xy_01 * x },
            { 2 * x, 4 * y }
        };
    }

    // добуток (A_0)^-1 * F_k
    static double[] FindMatrixMult(double[,] A, double[] F)
    {
        double detA = A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0];
        if (Math.Abs(detA) < 1e-10)
        {
            Console.WriteLine("The matrix is nearly singular, and the system cannot be solved.");
            Environment.Exit(0);
        }

        double[,] A_inv = new double[,] {
            { A[1, 1] / detA, -A[0, 1] / detA },
            { -A[1, 0] / detA, A[0, 0] / detA }
        };

        double[] result = new double[2];
        result[0] = A_inv[0, 0] * F[0] + A_inv[0, 1] * F[1];
        result[1] = A_inv[1, 0] * F[0] + A_inv[1, 1] * F[1];

        return result;
    }

    static void ModifiedNewtonMethod(double x, double y, double precision)
    {
        int maxIterations = 6;
        int k = 1;

        // ініціалізація x_0, A_0
        double[] xPrev = { x, y };
        double[,] A0 = JacobiMatrix(x, y);

        while (k < maxIterations)
        {
            // знаходження F_k
            double[] F = SystemOfEquations(x, y);

            // знаходження результату добутку (A_0)^-1 * F_k, що потрібно для ітераційного процесу
            double[] mult = FindMatrixMult(A0, F);

            // знаходження x_k+1 (x_k+1 = x_k - (A_0)^(-1) * F(x_k))
            x -= mult[0];
            y -= mult[1];

            // норма (максимальний за модулем елемент вектора)
            double norm = Math.Max(Math.Abs(x - xPrev[0]), Math.Abs(y - xPrev[1]));

            Console.WriteLine($"Iteration {k}:");
            Console.WriteLine($"F(x, y) = [{F[0]}, {F[1]}]");
            Console.WriteLine($"Norm = {norm}");
            Console.WriteLine($"\nx = {x}, y = {y}");
            Console.WriteLine("-----------------------");

            // перевірка умови припинення ітераційного процесу
            if (norm <= precision)
            {
                Console.WriteLine("Solution found:");
                Console.WriteLine($"x = {x}");
                Console.WriteLine($"y = {y}");
                Console.WriteLine($"Precision achieved: {norm}");
                return;
            }

            xPrev[0] = x;
            xPrev[1] = y;

            k++;
        }

        Console.WriteLine("No result found within five iterations.");
    }

    static void Main()
    {
        double x, y, precision;

        Console.WriteLine("System of equations:");
        Console.WriteLine("1. tan(xy + 0.1) = x^2");
        Console.WriteLine("2. x^2 + 2 * y^2 = 1");

        Console.WriteLine("\nEnter an initial approximation for vector x_0.");
        Console.Write("x: ");
        while (!double.TryParse(Console.ReadLine(), out x))
        {
            Console.WriteLine("Invalid input. Please enter a valid double number: ");
        }

        Console.Write("y: ");
        while (!double.TryParse(Console.ReadLine(), out y))
        {
            Console.WriteLine("Invalid input. Please enter a valid double number: ");
        }

        Console.Write("\nEnter a precision E: ");
        while (!double.TryParse(Console.ReadLine(), out precision))
        {
            Console.WriteLine("Invalid input. Please enter a valid number: ");
        }

        ModifiedNewtonMethod(x, y, precision);
    }
}