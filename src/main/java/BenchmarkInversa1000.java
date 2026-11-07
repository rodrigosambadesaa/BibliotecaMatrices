import ch.obermuhlner.math.big.BigComplex;
import comun.MatrizCompleja;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Locale;

/**
 * Reproducible, standalone benchmark for the library's exact inversa() API.
 * Run with: java -cp 'target/classes:target/dependency/*' BenchmarkInversa1000 array 1000
 */
public final class BenchmarkInversa1000 {
    private BenchmarkInversa1000() { }

    public static void main(String[] args) {
        if (args.length != 2) {
            System.err.println("Usage: BenchmarkInversa1000 <array|lista> <dimension>");
            System.exit(2);
        }
        String implementation = args[0];
        int n = Integer.parseInt(args[1]);
        if (n < 2 || n > 5000 || (!implementation.equals("array") && !implementation.equals("lista"))) {
            throw new IllegalArgumentException("Expected array|lista and 2 <= dimension <= 5000");
        }

        long setupStart = System.nanoTime();
        BigComplex[][] values = new BigComplex[n][n];
        BigComplex minusOne = BigComplex.valueOf(BigDecimal.ONE.negate());
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    values[i][j] = BigComplex.valueOf(BigDecimal.valueOf(4L * n + i % 5));
                } else {
                    int value = Math.floorMod(17 * i + 29 * j + 7, 3) - 1;
                    values[i][j] = value == 0 ? BigComplex.ZERO
                            : value == 1 ? BigComplex.ONE : minusOne;
                }
            }
        }
        MatrizCompleja original = implementation.equals("array")
                ? new array.Matriz(values) : new lista.Matriz(values);
        long setupNanos = System.nanoTime() - setupStart;
        values = null; // Avoid holding an unnecessary second reference to the input data.

        System.out.printf(Locale.ROOT, "BENCH_BEGIN impl=%s n=%d java=%s processors=%d max_heap_mb=%d%n",
                implementation, n, System.getProperty("java.version"),
                Runtime.getRuntime().availableProcessors(), Runtime.getRuntime().maxMemory() / 1048576);
        System.out.flush();

        long start = System.nanoTime();
        MatrizCompleja inverse = implementation.equals("array")
                ? ((array.Matriz) original).inversa() : ((lista.Matriz) original).inversa();
        long inverseNanos = System.nanoTime() - start;

        // Check nine selected cells of A*A^-1 and A^-1*A, without charging
        // an additional O(n^3) full-matrix multiplication to the measured time.
        MathContext mc = MathContext.DECIMAL128;
        BigDecimal maxError = BigDecimal.ZERO;
        int[] samples = {0, n / 2, n - 1};
        for (boolean reverse : new boolean[] {false, true}) {
            for (int i : samples) {
                for (int j : samples) {
                    BigComplex sum = BigComplex.ZERO;
                    for (int k = 0; k < n; k++) {
                        BigComplex left = reverse ? inverse.get(i, k) : original.get(i, k);
                        BigComplex right = reverse ? original.get(k, j) : inverse.get(k, j);
                        sum = sum.add(left.multiply(right, mc), mc);
                    }
                    BigComplex desired = i == j ? BigComplex.ONE : BigComplex.ZERO;
                    BigDecimal error = sum.subtract(desired, mc).abs(mc);
                    maxError = maxError.max(error);
                }
            }
        }
        System.out.printf(Locale.ROOT,
                "BENCH_RESULT impl=%s n=%d construct_ms=%.3f inverse_ms=%.3f sample_max_error=%s%n",
                implementation, n, setupNanos / 1_000_000.0,
                inverseNanos / 1_000_000.0, maxError.toEngineeringString());
        System.out.flush();
        if (maxError.compareTo(new BigDecimal("1E-18")) > 0) {
            throw new AssertionError("Inverse residual larger than 1E-18: " + maxError);
        }
    }
}
