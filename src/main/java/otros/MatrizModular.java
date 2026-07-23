package otros;

import comun.excepciones.DimensionMatrizException;
import comun.excepciones.MatrizSingularException;
import comun.excepciones.OperacionNoDefinidaException;

import java.util.Arrays;
import java.util.Objects;

/**
 * Matriz sobre el cuerpo primo F_p. El constructor verifica primalidad y todas
 * las operaciones normalizan residuos a [0,p).
 */
public final class MatrizModular {
    private final int[][] datos;
    private final int modulo;

    public MatrizModular(int[][] entrada, int modulo) {
        Objects.requireNonNull(entrada);
        if (!esPrimo(modulo)) throw new OperacionNoDefinidaException("El módulo debe ser primo");
        if (entrada.length == 0 || entrada[0] == null || entrada[0].length == 0)
            throw new DimensionMatrizException("Matriz modular vacía");
        this.modulo = modulo;
        datos = new int[entrada.length][entrada[0].length];
        for (int i = 0; i < entrada.length; i++) {
            if (entrada[i] == null || entrada[i].length != datos[0].length)
                throw new DimensionMatrizException("Matriz modular irregular");
            for (int j = 0; j < datos[i].length; j++) datos[i][j] = normalizar(entrada[i][j]);
        }
    }
    public int filas() { return datos.length; }
    public int columnas() { return datos[0].length; }
    public int modulo() { return modulo; }
    public int get(int i, int j) { return datos[i][j]; }
    public int[][] getMatriz() {
        int[][] c = new int[filas()][];
        for (int i = 0; i < filas(); i++) c[i] = datos[i].clone();
        return c;
    }
    public MatrizModular suma(MatrizModular otra) {
        exigirCompatible(otra, true); int[][] r = new int[filas()][columnas()];
        for (int i = 0; i < filas(); i++) for (int j = 0; j < columnas(); j++) r[i][j] = normalizar(datos[i][j] + otra.datos[i][j]);
        return new MatrizModular(r, modulo);
    }
    public MatrizModular producto(MatrizModular otra) {
        exigirCompatible(otra, false);
        if (columnas() != otra.filas()) throw new DimensionMatrizException("Producto modular incompatible");
        int[][] r = new int[filas()][otra.columnas()];
        for (int i = 0; i < filas(); i++) for (int j = 0; j < otra.columnas(); j++)
            for (int k = 0; k < columnas(); k++) r[i][j] = normalizar(r[i][j] + datos[i][k] * otra.datos[k][j]);
        return new MatrizModular(r, modulo);
    }
    public int determinante() {
        if (filas() != columnas()) throw new DimensionMatrizException("Determinante modular exige matriz cuadrada");
        int[][] a = getMatriz(); int det = 1;
        for (int k = 0; k < filas(); k++) {
            int p = k; while (p < filas() && a[p][k] == 0) p++;
            if (p == filas()) return 0;
            if (p != k) { int[] t = a[p]; a[p] = a[k]; a[k] = t; det = normalizar(-det); }
            det = normalizar(det * a[k][k]);
            int inv = inverso(a[k][k]);
            for (int i = k + 1; i < filas(); i++) {
                int f = normalizar(a[i][k] * inv);
                for (int j = k; j < columnas(); j++) a[i][j] = normalizar(a[i][j] - f * a[k][j]);
            }
        }
        return det;
    }
    public int rango() {
        return gauss().pivotes;
    }
    public MatrizModular gaussReducida() { return new MatrizModular(gauss().datos, modulo); }
    public MatrizModular inversa() {
        if (filas() != columnas()) throw new DimensionMatrizException("Inversa modular exige matriz cuadrada");
        int n = filas(); int[][] aug = new int[n][2 * n];
        for (int i = 0; i < n; i++) { System.arraycopy(datos[i], 0, aug[i], 0, n); aug[i][n + i] = 1; }
        ResultadoGauss g = gauss(aug);
        if (g.pivotes != n) throw new MatrizSingularException("Matriz singular módulo " + modulo);
        int[][] inv = new int[n][n];
        for (int i = 0; i < n; i++) System.arraycopy(g.datos[i], n, inv[i], 0, n);
        return new MatrizModular(inv, modulo);
    }
    private ResultadoGauss gauss() { return gauss(getMatriz()); }
    private ResultadoGauss gauss(int[][] a) {
        int fila = 0;
        for (int col = 0; col < a[0].length && fila < a.length; col++) {
            int p = fila; while (p < a.length && a[p][col] == 0) p++;
            if (p == a.length) continue;
            int[] t = a[p]; a[p] = a[fila]; a[fila] = t;
            int inv = inverso(a[fila][col]);
            for (int j = col; j < a[0].length; j++) a[fila][j] = normalizar(a[fila][j] * inv);
            for (int i = 0; i < a.length; i++) if (i != fila) {
                int f = a[i][col];
                for (int j = col; j < a[0].length; j++) a[i][j] = normalizar(a[i][j] - f * a[fila][j]);
            }
            fila++;
        }
        return new ResultadoGauss(a, fila);
    }
    private int inverso(int a) {
        int t = 0, nt = 1, r = modulo, nr = normalizar(a);
        while (nr != 0) { int q = r / nr; int x = t - q * nt; t = nt; nt = x; x = r - q * nr; r = nr; nr = x; }
        if (r != 1) throw new OperacionNoDefinidaException("Elemento no invertible");
        return normalizar(t);
    }
    private int normalizar(int x) { return Math.floorMod(x, modulo); }
    private void exigirCompatible(MatrizModular o, boolean dimensiones) {
        if (o == null || modulo != o.modulo || (dimensiones && (filas() != o.filas() || columnas() != o.columnas())))
            throw new DimensionMatrizException("Matrices modulares incompatibles");
    }
    private static boolean esPrimo(int n) {
        if (n < 2) return false;
        for (int i = 2; (long) i * i <= n; i++) if (n % i == 0) return false;
        return true;
    }
    private record ResultadoGauss(int[][] datos, int pivotes) {}
    @Override public boolean equals(Object o) {
        return o instanceof MatrizModular m && modulo == m.modulo && Arrays.deepEquals(datos, m.datos);
    }
    @Override public int hashCode() { return 31 * modulo + Arrays.deepHashCode(datos); }
    @Override public String toString() { return Arrays.deepToString(datos) + " (mod " + modulo + ')'; }
}
