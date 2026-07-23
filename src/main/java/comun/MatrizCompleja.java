package comun;

import ch.obermuhlner.math.big.BigComplex;

import java.math.BigDecimal;
import java.math.MathContext;

/** Vista pública de solo lectura que permite compartir algoritmos sin compartir almacenamiento. */
public interface MatrizCompleja {
    int filas();
    int columnas();
    BigComplex get(int fila, int columna);
    MathContext mathContext();
    BigDecimal tolerancia();
    array.Matriz comoMatrizArray();
    lista.Matriz comoMatrizLista();
}
