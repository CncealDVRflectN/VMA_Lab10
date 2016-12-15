package by.vma.lab10;

public class Main {
    private static class Matrix {
        public double[][] matrix;
        private int lines;
        private int columns;

        public Matrix(int lines, int columns) throws Exception {
            if (lines < 1 || columns < 1) {
                throw new Exception("Неверный размер.");
            }
            this.lines = lines;
            this.columns = columns;
            this.matrix = new double[lines][columns];
        }

        public Matrix(Matrix init) throws Exception {
            this(init.getLines(), init.getColumns());
            for (int i = 0; i < lines; i++) {
                for (int j = 0; j < columns; j++) {
                    this.matrix[i][j] = init.matrix[i][j];
                }
            }
        }

        public int getLines() {
            return lines;
        }

        public int getColumns() {
            return columns;
        }

        public Vector getColumn(int index) throws Exception {
            if(index < 0 || index >= columns) {
                throw new Exception("Неверный индекс.");
            }
            Vector result = new Vector(n);
            for(int i = 0; i < n; i++) {
                result.vector[i] = this.matrix[i][index];
            }
            return result;
        }

        public void print() {
            for (double[] i : matrix) {
                for (double j : i) {
                    System.out.printf("%.5f", j);
                    System.out.print("  ");
                }
                System.out.println();
            }
        }

        public void swap(int fi, int fj, int si, int sj) {
            double tmp = matrix[fi][fj];
            matrix[fi][fj] = matrix[si][sj];
            matrix[si][sj] = tmp;
        }

        public void fillDefault() {
            double[][] a = {{0.6444, 0.0000, -0.1683, 0.1184, 0.1973},
                    {-0.0395, 0.4208, 0.0000, -0.0802, 0.0263},
                    {0.0132, -0.1184, 0.7627, 0.0145, 0.0460},
                    {0.0395, 0.0000, -0.0960, 0.7627, 0.0000},
                    {0.0263, -0.0395, 0.1907, -0.0158, 0.5523}};
            this.lines = 5;
            this.columns = 5;
            this.matrix = a;
        }

        public void fillE(int n) {
            lines = n;
            columns = n;
            matrix = new double[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (i == j) {
                        matrix[i][j] = 1;
                    } else {
                        matrix[i][j] = 0;
                    }
                }
            }
        }

        public Vector mul(Vector vector) throws Exception {
            if (columns != vector.getLength()) {
                throw new Exception("Неверная матрица или вектор.");
            }
            Vector result = new Vector(vector.getLength());
            for (int i = 0; i < lines; i++) {
                result.vector[i] = 0;
                for (int j = 0; j < columns; j++) {
                    result.vector[i] += matrix[i][j] * vector.vector[j];
                }
            }
            return result;
        }

        public Matrix mul(Matrix mtr) throws Exception {
            if (columns != mtr.getLines()) {
                throw new Exception("Неверная матрица.");
            }
            Matrix result = new Matrix(lines, mtr.getColumns());
            for (int i = 0; i < result.getLines(); i++) {
                for (int j = 0; j < result.getColumns(); j++) {
                    result.matrix[i][j] = 0;
                    for (int k = 0; k < columns; k++) {
                        result.matrix[i][j] += this.matrix[i][k] * mtr.matrix[k][j];
                    }
                }
            }
            return result;
        }

        public Matrix transpose() throws Exception {
            if (lines != columns) {
                throw new Exception("Неверная матрица.");
            }
            Matrix result = new Matrix(this);
            for (int i = 0; i < lines; i++) {
                for (int j = i + 1; j < columns; j++) {
                    result.swap(i, j, j, i);
                }
            }
            return result;
        }
    }

    private static class Vector {
        public double[] vector;
        private int length;

        public Vector(int length) throws Exception {
            if (length < 1) {
                throw new Exception("Неверный размер.");
            }
            this.length = length;
            vector = new double[length];
        }

        public int getLength() {
            return length;
        }

        public void print(boolean exponent) {
            for (double item : vector) {
                if (exponent) {
                    System.out.printf("%e\n", item);
                } else {
                    System.out.printf("%.5f\n", item);
                }
            }
        }

        public Vector mul(double num) throws Exception {
            Vector result = new Vector(length);
            for (int i = 0; i < length; i++) {
                result.vector[i] = this.vector[i] * num;
            }
            return result;
        }

        public Vector subtract(Vector sub) throws Exception {
            if (length != sub.getLength()) {
                throw new Exception("Неверный вектор.");
            }
            Vector result = new Vector(length);
            for (int i = 0; i < length; i++) {
                result.vector[i] = this.vector[i] - sub.vector[i];
            }
            return result;
        }

        public double normI() {
            double max = Math.abs(vector[0]);
            for (int i = 1; i < length; i++) {
                if (Math.abs(vector[i]) > max) {
                    max = Math.abs(vector[i]);
                }
            }
            return max;
        }
    }

    private static Matrix A;
    private static Matrix U;
    private static final int n = 5;
    private static final double epsilon = 0.00001;

    public static void main(String[] args) {
        Vector result;
        Vector eigen;
        Vector r;
        try {
            A = new Matrix(n, n);
            U = new Matrix(n, n);
            A.fillDefault();
            A = A.transpose().mul(A);
            System.out.println("Матрица A:");
            A.print();
            System.out.println();
            U.fillE(n);
            System.out.println("k: " + ((Math.log(epsilon) - Math.log(calculateT())) / (Math.log(1 - 2.0 / (n * (n - 1))))));
            result = rotationMethod();
            A.fillDefault();
            A = A.transpose().mul(A);
            for(int i = 0; i < n; i++) {
                System.out.println("Собственное значение " + (i + 1) + ": " + result.vector[i]);
                eigen = U.getColumn(i);
                System.out.println("Соответствующий собственный вектор:");
                eigen.print(false);
                r = A.mul(eigen).subtract(eigen.mul(result.vector[i]));
                System.out.println("Вектор невязки:");
                r.print(true);
                System.out.println("Норма вектора невязки: " + r.normI());
                System.out.println();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static Vector rotationMethod() throws Exception {
        Vector result = new Vector(n);
        double cos;
        double sin;
        double cosd;
        double tgd;
        double prevA;
        double prevE;
        double max;
        int maxl;
        int maxm;
        int counter = 0;
        do {
            max = A.matrix[0][1];
            maxl = 0;
            maxm = 1;
            for (int l = 0; l < n; l++) {
                for (int m = l + 1; m < n; m++) {
                    if (Math.abs(A.matrix[l][m]) > Math.abs(max)) {
                        max = A.matrix[l][m];
                        maxl = l;
                        maxm = m;
                    }
                }
            }
            tgd = 2 * max / (A.matrix[maxl][maxl] - A.matrix[maxm][maxm]);
            cosd = 1 / Math.sqrt(1 + Math.pow(tgd, 2));
            cos = Math.sqrt((1 + cosd) / 2);
            sin = Math.sqrt((1 - cosd) / 2) * Math.signum(tgd);
            for (int i = 0; i < n; i++) {
                prevA = A.matrix[i][maxl];
                prevE = U.matrix[i][maxl];
                A.matrix[i][maxl] = prevA * cos + A.matrix[i][maxm] * sin;
                A.matrix[i][maxm] = -prevA * sin + A.matrix[i][maxm] * cos;
                U.matrix[i][maxl] = prevE * cos + U.matrix[i][maxm] * sin;
                U.matrix[i][maxm] = -prevE * sin + U.matrix[i][maxm] * cos;
            }
            for (int i = 0; i < n; i++) {
                prevA = A.matrix[maxl][i];
                A.matrix[maxl][i] = prevA * cos + A.matrix[maxm][i] * sin;
                A.matrix[maxm][i] = -prevA * sin + A.matrix[maxm][i] * cos;
            }
            counter++;
        } while (calculateT() > epsilon);
        for(int i = 0; i < n; i++){
            result.vector[i] = A.matrix[i][i];
        }
        System.out.println("Количество итераций: " + counter);
        System.out.println();
        return result;
    }

    private static double calculateT() {
        double result = 0;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (i != j) {
                    result += Math.pow(A.matrix[i][j], 2);
                }
            }
        }
        return 2 * result;
    }
}
