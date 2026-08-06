package comun;

import ch.obermuhlner.math.big.BigComplex;

import java.math.BigDecimal;
import java.math.MathContext;

/**
 * Vista de solo lectura de una matriz compleja tridimensional.
 * El orden de los ejes es siempre {@code [capa][fila][columna]}.
 */
public interface Matriz3DCompleja {
    int capas();
    int filas();
    int columnas();
    BigComplex get(int capa, int fila, int columna);
    MathContext mathContext();
    BigDecimal tolerancia();
    array.Matriz3D comoMatriz3DArray();
    lista.Matriz3D comoMatriz3DLista();
}
