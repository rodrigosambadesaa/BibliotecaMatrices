package otros;

import comun.excepciones.DimensionMatrizException;

import java.util.Arrays;
import java.util.Objects;

/** Matriz sobre el semianillo booleano; no reutiliza aritmética compleja. */
public final class MatrizBooleana {
    private final boolean[][] datos;

    public MatrizBooleana(boolean[][] entrada) {
        Objects.requireNonNull(entrada);
        if (entrada.length == 0 || entrada[0] == null || entrada[0].length == 0)
            throw new DimensionMatrizException("Matriz booleana vacía");
        datos = new boolean[entrada.length][entrada[0].length];
        for (int i = 0; i < entrada.length; i++) {
            if (entrada[i] == null || entrada[i].length != datos[0].length)
                throw new DimensionMatrizException("Matriz booleana irregular");
            datos[i] = entrada[i].clone();
        }
    }
    public int filas() { return datos.length; }
    public int columnas() { return datos[0].length; }
    public boolean get(int i, int j) { return datos[i][j]; }
    public boolean[][] getMatriz() {
        boolean[][] c = new boolean[filas()][];
        for (int i = 0; i < filas(); i++) c[i] = datos[i].clone();
        return c;
    }
    public MatrizBooleana suma(MatrizBooleana otra) {
        exigirIgual(otra); boolean[][] r = new boolean[filas()][columnas()];
        for (int i = 0; i < filas(); i++) for (int j = 0; j < columnas(); j++) r[i][j] = datos[i][j] || otra.datos[i][j];
        return new MatrizBooleana(r);
    }
    public MatrizBooleana producto(MatrizBooleana otra) {
        if (columnas() != otra.filas()) throw new DimensionMatrizException("Producto booleano incompatible");
        boolean[][] r = new boolean[filas()][otra.columnas()];
        for (int i = 0; i < filas(); i++) for (int j = 0; j < otra.columnas(); j++)
            for (int k = 0; k < columnas(); k++) r[i][j] |= datos[i][k] && otra.datos[k][j];
        return new MatrizBooleana(r);
    }
    public MatrizBooleana clausuraTransitiva() {
        if (filas() != columnas()) throw new DimensionMatrizException("La clausura exige matriz cuadrada");
        boolean[][] r = getMatriz();
        for (int i = 0; i < filas(); i++) r[i][i] = true;
        for (int k = 0; k < filas(); k++) for (int i = 0; i < filas(); i++) for (int j = 0; j < filas(); j++)
            r[i][j] |= r[i][k] && r[k][j];
        return new MatrizBooleana(r);
    }
    private void exigirIgual(MatrizBooleana o) {
        if (o == null || filas() != o.filas() || columnas() != o.columnas()) throw new DimensionMatrizException("Dimensiones distintas");
    }
    @Override public boolean equals(Object o) { return o instanceof MatrizBooleana m && Arrays.deepEquals(datos, m.datos); }
    @Override public int hashCode() { return Arrays.deepHashCode(datos); }
    @Override public String toString() { return Arrays.deepToString(datos); }
}
