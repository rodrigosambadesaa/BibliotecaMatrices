package comun;

import ch.obermuhlner.math.big.BigComplex;
import comun.excepciones.DimensionMatrizException;
import comun.excepciones.DivisionPorCeroException;
import comun.excepciones.IndiceMatrizException;
import comun.excepciones.OperacionNoDefinidaException;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.function.BiFunction;
import java.util.function.BinaryOperator;
import java.util.function.UnaryOperator;

/**
 * Núcleo inmutable común para matrices 3D complejas. Las operaciones que
 * tienen significado matricial se ejecutan independientemente por capa; las
 * demás actúan sobre los tres ejes completos.
 *
 * @param <T> implementación tridimensional concreta
 * @param <P> implementación bidimensional usada para los cortes
 */
public abstract class AbstractMatriz3D<T extends AbstractMatriz3D<T, P>, P extends AbstractMatriz<P>>
        implements Matriz3DCompleja {
    private final int capas;
    private final int filas;
    private final int columnas;
    private final MathContext mathContext;
    private final BigDecimal tolerancia;

    protected AbstractMatriz3D(int capas, int filas, int columnas,
                               MathContext mathContext, BigDecimal tolerancia) {
        if (capas <= 0 || filas <= 0 || columnas <= 0) {
            throw new DimensionMatrizException("Las dimensiones 3D deben ser positivas: "
                    + capas + "x" + filas + "x" + columnas);
        }
        this.capas = capas;
        this.filas = filas;
        this.columnas = columnas;
        this.mathContext = Objects.requireNonNull(mathContext, "El MathContext no puede ser null");
        Complejos.validarTolerancia(tolerancia);
        this.tolerancia = tolerancia;
    }

    protected abstract T crear(BigComplex[][][] datos, MathContext contexto, BigDecimal tolerancia);
    protected abstract P crearPlano(BigComplex[][] datos, MathContext contexto, BigDecimal tolerancia);

    @Override public final int capas() { return capas; }
    @Override public final int filas() { return filas; }
    @Override public final int columnas() { return columnas; }
    public final int getCapas() { return capas; }
    public final int getFilas() { return filas; }
    public final int getColumnas() { return columnas; }
    public final int numeroElementos() { return Math.multiplyExact(capas, Math.multiplyExact(filas, columnas)); }
    @Override public final MathContext mathContext() { return mathContext; }
    @Override public final BigDecimal tolerancia() { return tolerancia; }
    public final BigComplex obtener(int capa, int fila, int columna) { return get(capa, fila, columna); }

    protected final void validarIndice(int capa, int fila, int columna) {
        if (capa < 0 || capa >= capas || fila < 0 || fila >= filas || columna < 0 || columna >= columnas) {
            throw new IndiceMatrizException("Índice (" + capa + ',' + fila + ',' + columna + ") fuera de "
                    + capas + "x" + filas + "x" + columnas);
        }
    }

    protected final BigComplex redondear(BigComplex valor) {
        return Objects.requireNonNull(valor, "El valor no puede ser null").round(mathContext);
    }

    protected final boolean esCero(BigComplex valor) {
        return Complejos.esCero(valor, tolerancia, mathContext);
    }

    protected final void mismaPolitica(Matriz3DCompleja otra) {
        Objects.requireNonNull(otra, "La matriz 3D no puede ser null");
        if (!mathContext.equals(otra.mathContext()) || tolerancia.compareTo(otra.tolerancia()) != 0) {
            throw new OperacionNoDefinidaException("No se mezclan políticas de precisión distintas");
        }
    }

    protected final void mismasDimensiones(Matriz3DCompleja otra) {
        mismaPolitica(otra);
        if (capas != otra.capas() || filas != otra.filas() || columnas != otra.columnas()) {
            throw new DimensionMatrizException("Dimensiones 3D incompatibles: " + dimensiones()
                    + " y " + otra.capas() + "x" + otra.filas() + "x" + otra.columnas());
        }
    }

    public final String dimensiones() { return capas + "x" + filas + "x" + columnas; }

    protected final BigComplex[][][] copiaDatos() {
        BigComplex[][][] copia = new BigComplex[capas][filas][columnas];
        for (int k = 0; k < capas; k++)
            for (int i = 0; i < filas; i++)
                for (int j = 0; j < columnas; j++) copia[k][i][j] = get(k, i, j);
        return copia;
    }

    public final BigComplex[][][] getMatriz() { return copiaDatos(); }

    @Override public final array.Matriz3D comoMatriz3DArray() {
        return new array.Matriz3D(this, mathContext, tolerancia);
    }

    @Override public final lista.Matriz3D comoMatriz3DLista() {
        return new lista.Matriz3D(this, mathContext, tolerancia);
    }

    /** Corte perpendicular al eje de capas, de tamaño filas × columnas. */
    public final P capa(int capa) {
        validarIndice(capa, 0, 0);
        BigComplex[][] datos = new BigComplex[filas][columnas];
        for (int i = 0; i < filas; i++)
            for (int j = 0; j < columnas; j++) datos[i][j] = get(capa, i, j);
        return crearPlano(datos, mathContext, tolerancia);
    }

    /** Corte perpendicular al eje de filas, de tamaño capas × columnas. */
    public final P cortePorFila(int fila) {
        validarIndice(0, fila, 0);
        BigComplex[][] datos = new BigComplex[capas][columnas];
        for (int k = 0; k < capas; k++)
            for (int j = 0; j < columnas; j++) datos[k][j] = get(k, fila, j);
        return crearPlano(datos, mathContext, tolerancia);
    }

    /** Corte perpendicular al eje de columnas, de tamaño capas × filas. */
    public final P cortePorColumna(int columna) {
        validarIndice(0, 0, columna);
        BigComplex[][] datos = new BigComplex[capas][filas];
        for (int k = 0; k < capas; k++)
            for (int i = 0; i < filas; i++) datos[k][i] = get(k, i, columna);
        return crearPlano(datos, mathContext, tolerancia);
    }

    /**
     * Aplica a cada capa cualquier operación 2D que devuelva otra matriz.
     * Todas las salidas deben compartir dimensiones y política de precisión.
     * Permite reutilizar, por ejemplo, inversas, RREF, exponenciales o
     * factorizaciones que se materialicen como matriz.
     */
    public final T aplicarPorCapas(UnaryOperator<P> operacion) {
        Objects.requireNonNull(operacion, "La operación no puede ser null");
        List<P> resultados = new ArrayList<>(capas);
        for (int k = 0; k < capas; k++) resultados.add(Objects.requireNonNull(operacion.apply(capa(k)),
                "La operación no puede devolver null"));
        return ensamblarResultados(resultados);
    }

    /**
     * Combina pares de capas con cualquier operación binaria 2D que devuelva
     * una matriz, como producto, resolución de sistemas o suma por bloques.
     */
    public final T combinarPorCapas(Matriz3DCompleja otra, BiFunction<P, P, P> operacion) {
        mismaPolitica(otra);
        if (capas != otra.capas()) throw new DimensionMatrizException("El número de capas debe coincidir");
        Objects.requireNonNull(operacion, "La operación no puede ser null");
        List<P> resultados = new ArrayList<>(capas);
        for (int k = 0; k < capas; k++) resultados.add(Objects.requireNonNull(
                operacion.apply(capa(k), capaDe(otra, k)), "La operación no puede devolver null"));
        return ensamblarResultados(resultados);
    }

    public final T conElemento(int capa, int fila, int columna, BigComplex valor) {
        validarIndice(capa, fila, columna);
        BigComplex[][][] datos = copiaDatos();
        datos[capa][fila][columna] = redondear(valor);
        return crear(datos, mathContext, tolerancia);
    }

    public final T conCapa(int indice, MatrizCompleja nuevaCapa) {
        validarIndice(indice, 0, 0);
        Objects.requireNonNull(nuevaCapa, "La capa no puede ser null");
        if (!mathContext.equals(nuevaCapa.mathContext()) || tolerancia.compareTo(nuevaCapa.tolerancia()) != 0)
            throw new OperacionNoDefinidaException("No se mezclan políticas de precisión distintas");
        if (nuevaCapa.filas() != filas || nuevaCapa.columnas() != columnas)
            throw new DimensionMatrizException("La capa debe medir " + filas + "x" + columnas);
        BigComplex[][][] datos = copiaDatos();
        for (int i = 0; i < filas; i++)
            for (int j = 0; j < columnas; j++) datos[indice][i][j] = nuevaCapa.get(i, j);
        return crear(datos, mathContext, tolerancia);
    }

    /** Submatriz con los seis límites desde inclusivos y hasta exclusivos. */
    public final T submatriz(int capaDesde, int capaHasta, int filaDesde, int filaHasta,
                             int columnaDesde, int columnaHasta) {
        if (capaDesde < 0 || capaHasta > capas || capaDesde >= capaHasta
                || filaDesde < 0 || filaHasta > filas || filaDesde >= filaHasta
                || columnaDesde < 0 || columnaHasta > columnas || columnaDesde >= columnaHasta) {
            throw new IndiceMatrizException("Intervalo de submatriz 3D inválido");
        }
        BigComplex[][][] datos = new BigComplex[capaHasta - capaDesde][filaHasta - filaDesde][columnaHasta - columnaDesde];
        for (int k = capaDesde; k < capaHasta; k++)
            for (int i = filaDesde; i < filaHasta; i++)
                for (int j = columnaDesde; j < columnaHasta; j++)
                    datos[k - capaDesde][i - filaDesde][j - columnaDesde] = get(k, i, j);
        return crear(datos, mathContext, tolerancia);
    }

    public final T extraerCapas(int... indices) {
        Objects.requireNonNull(indices, "Los índices no pueden ser null");
        if (indices.length == 0) throw new DimensionMatrizException("Debe extraerse alguna capa");
        BigComplex[][][] datos = new BigComplex[indices.length][filas][columnas];
        for (int n = 0; n < indices.length; n++) {
            validarIndice(indices[n], 0, 0);
            for (int i = 0; i < filas; i++)
                for (int j = 0; j < columnas; j++) datos[n][i][j] = get(indices[n], i, j);
        }
        return crear(datos, mathContext, tolerancia);
    }

    public final T insertarCapa(int posicion, MatrizCompleja nuevaCapa) {
        if (posicion < 0 || posicion > capas) throw new IndiceMatrizException("Posición de capa inválida");
        Objects.requireNonNull(nuevaCapa, "La capa no puede ser null");
        if (!mathContext.equals(nuevaCapa.mathContext()) || tolerancia.compareTo(nuevaCapa.tolerancia()) != 0)
            throw new OperacionNoDefinidaException("No se mezclan políticas de precisión distintas");
        if (nuevaCapa.filas() != filas || nuevaCapa.columnas() != columnas)
            throw new DimensionMatrizException("La capa debe medir " + filas + "x" + columnas);
        BigComplex[][][] datos = new BigComplex[capas + 1][filas][columnas];
        for (int k = 0, origen = 0; k < datos.length; k++) {
            if (k == posicion) {
                for (int i = 0; i < filas; i++)
                    for (int j = 0; j < columnas; j++) datos[k][i][j] = nuevaCapa.get(i, j);
            } else {
                for (int i = 0; i < filas; i++)
                    for (int j = 0; j < columnas; j++) datos[k][i][j] = get(origen, i, j);
                origen++;
            }
        }
        return crear(datos, mathContext, tolerancia);
    }

    public final T eliminarCapa(int indice) {
        if (capas == 1) throw new DimensionMatrizException("No se admite una matriz 3D sin capas");
        validarIndice(indice, 0, 0);
        int[] indices = new int[capas - 1];
        for (int k = 0, n = 0; k < capas; k++) if (k != indice) indices[n++] = k;
        return extraerCapas(indices);
    }

    /** Concatena por el eje indicado: 0 capas, 1 filas o 2 columnas. */
    public final T concatenar(Matriz3DCompleja otra, int eje) {
        mismaPolitica(otra);
        validarEje(eje);
        int[] a = {capas, filas, columnas};
        int[] b = {otra.capas(), otra.filas(), otra.columnas()};
        for (int x = 0; x < 3; x++)
            if (x != eje && a[x] != b[x]) throw new DimensionMatrizException("Dimensiones incompatibles para concatenar");
        int[] salida = a.clone();
        salida[eje] = Math.addExact(a[eje], b[eje]);
        BigComplex[][][] datos = new BigComplex[salida[0]][salida[1]][salida[2]];
        for (int k = 0; k < salida[0]; k++) for (int i = 0; i < salida[1]; i++) for (int j = 0; j < salida[2]; j++) {
            int[] p = {k, i, j};
            boolean segundo = p[eje] >= a[eje];
            if (segundo) p[eje] -= a[eje];
            datos[k][i][j] = segundo ? otra.get(p[0], p[1], p[2]) : get(p[0], p[1], p[2]);
        }
        return crear(datos, mathContext, tolerancia);
    }

    public final T concatenarCapas(Matriz3DCompleja otra) { return concatenar(otra, 0); }
    public final T concatenarFilas(Matriz3DCompleja otra) { return concatenar(otra, 1); }
    public final T concatenarColumnas(Matriz3DCompleja otra) { return concatenar(otra, 2); }

    public final T redimensionar(int nuevasCapas, int nuevasFilas, int nuevasColumnas, BigComplex relleno) {
        if (nuevasCapas <= 0 || nuevasFilas <= 0 || nuevasColumnas <= 0)
            throw new DimensionMatrizException("Las dimensiones 3D deben ser positivas");
        BigComplex valor = redondear(relleno);
        BigComplex[][][] datos = nuevosDatos(nuevasCapas, nuevasFilas, nuevasColumnas, valor);
        for (int k = 0; k < Math.min(capas, nuevasCapas); k++)
            for (int i = 0; i < Math.min(filas, nuevasFilas); i++)
                for (int j = 0; j < Math.min(columnas, nuevasColumnas); j++) datos[k][i][j] = get(k, i, j);
        return crear(datos, mathContext, tolerancia);
    }

    public final T rellenar(BigComplex valor) {
        return crear(nuevosDatos(capas, filas, columnas, redondear(valor)), mathContext, tolerancia);
    }

    public final T map(UnaryOperator<BigComplex> funcion) {
        Objects.requireNonNull(funcion, "La función no puede ser null");
        BigComplex[][][] datos = new BigComplex[capas][filas][columnas];
        for (int k = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            datos[k][i][j] = redondear(funcion.apply(get(k, i, j)));
        return crear(datos, mathContext, tolerancia);
    }

    public final T zip(Matriz3DCompleja otra, BinaryOperator<BigComplex> funcion) {
        mismasDimensiones(otra);
        Objects.requireNonNull(funcion, "La función no puede ser null");
        BigComplex[][][] datos = new BigComplex[capas][filas][columnas];
        for (int k = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            datos[k][i][j] = redondear(funcion.apply(get(k, i, j), otra.get(k, i, j)));
        return crear(datos, mathContext, tolerancia);
    }

    public final T suma(Matriz3DCompleja otra) { return zip(otra, (a, b) -> a.add(b, mathContext)); }
    public final T resta(Matriz3DCompleja otra) { return zip(otra, (a, b) -> a.subtract(b, mathContext)); }
    public final T opuesta() { return map(BigComplex::negate); }
    public final T conjugada() { return map(BigComplex::conjugate); }
    public final T productoEscalar(BigComplex escalar) {
        Objects.requireNonNull(escalar, "El escalar no puede ser null");
        return map(x -> x.multiply(escalar, mathContext));
    }
    public final T productoEscalar(BigDecimal escalar) { return productoEscalar(BigComplex.valueOf(Objects.requireNonNull(escalar))); }
    public final T productoEscalar(double escalar) { return productoEscalar(Complejos.de(escalar)); }
    public final T dividirPorEscalar(BigComplex escalar) {
        Objects.requireNonNull(escalar, "El escalar no puede ser null");
        if (esCero(escalar)) throw new DivisionPorCeroException("No se puede dividir por cero");
        return map(x -> x.divide(escalar, mathContext));
    }
    public final T dividirPorEscalar(BigDecimal escalar) { return dividirPorEscalar(BigComplex.valueOf(Objects.requireNonNull(escalar))); }
    public final T productoHadamard(Matriz3DCompleja otra) { return zip(otra, (a, b) -> a.multiply(b, mathContext)); }
    public final T cocienteHadamard(Matriz3DCompleja otra) {
        mismasDimensiones(otra);
        for (int k = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            if (esCero(otra.get(k, i, j))) throw new DivisionPorCeroException("Divisor Hadamard nulo");
        return zip(otra, (a, b) -> a.divide(b, mathContext));
    }

    /** Multiplicación matricial independiente para cada capa. */
    public final T productoPorCapas(Matriz3DCompleja otra) {
        mismaPolitica(otra);
        if (capas != otra.capas() || columnas != otra.filas())
            throw new DimensionMatrizException("Producto por capas incompatible");
        BigComplex[][][] datos = nuevosDatos(capas, filas, otra.columnas(), BigComplex.ZERO);
        for (int k = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < otra.columnas(); j++) {
            BigComplex suma = BigComplex.ZERO;
            for (int n = 0; n < columnas; n++)
                suma = suma.add(get(k, i, n).multiply(otra.get(k, n, j), mathContext), mathContext);
            datos[k][i][j] = suma.round(mathContext);
        }
        return crear(datos, mathContext, tolerancia);
    }

    /** Alias del producto matricial por capas. */
    public final T producto(Matriz3DCompleja otra) { return productoPorCapas(otra); }

    /** Multiplica todas las capas por la misma matriz bidimensional derecha. */
    public final T productoPorCapas(MatrizCompleja derecha) {
        Objects.requireNonNull(derecha, "La matriz no puede ser null");
        if (!mathContext.equals(derecha.mathContext()) || tolerancia.compareTo(derecha.tolerancia()) != 0)
            throw new OperacionNoDefinidaException("No se mezclan políticas de precisión distintas");
        if (columnas != derecha.filas()) throw new DimensionMatrizException("Producto por capas incompatible");
        BigComplex[][][] datos = nuevosDatos(capas, filas, derecha.columnas(), BigComplex.ZERO);
        for (int k = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < derecha.columnas(); j++) {
            BigComplex suma = BigComplex.ZERO;
            for (int n = 0; n < columnas; n++)
                suma = suma.add(get(k, i, n).multiply(derecha.get(n, j), mathContext), mathContext);
            datos[k][i][j] = suma.round(mathContext);
        }
        return crear(datos, mathContext, tolerancia);
    }

    /** Transpone filas y columnas dentro de cada capa. */
    public final T traspuestaPorCapas() {
        BigComplex[][][] datos = new BigComplex[capas][columnas][filas];
        for (int k = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            datos[k][j][i] = get(k, i, j);
        return crear(datos, mathContext, tolerancia);
    }

    public final T adjuntaHermiticaPorCapas() { return conjugada().traspuestaPorCapas(); }
    public final T traspuesta() { return traspuestaPorCapas(); }
    public final T adjuntaHermitica() { return adjuntaHermiticaPorCapas(); }

    /**
     * Permuta los ejes. Cada argumento indica el eje de entrada que ocupará
     * esa posición de salida; por ejemplo {@code (1, 0, 2)} intercambia capas y filas.
     */
    public final T permutarEjes(int eje0, int eje1, int eje2) {
        int[] permutacion = {eje0, eje1, eje2};
        validarPermutacion(permutacion);
        int[] entrada = {capas, filas, columnas};
        int[] salida = {entrada[eje0], entrada[eje1], entrada[eje2]};
        BigComplex[][][] datos = new BigComplex[salida[0]][salida[1]][salida[2]];
        for (int a = 0; a < salida[0]; a++) for (int b = 0; b < salida[1]; b++) for (int c = 0; c < salida[2]; c++) {
            int[] origen = new int[3];
            origen[eje0] = a; origen[eje1] = b; origen[eje2] = c;
            datos[a][b][c] = get(origen[0], origen[1], origen[2]);
        }
        return crear(datos, mathContext, tolerancia);
    }

    public final BigComplex sumaElementos() {
        BigComplex suma = BigComplex.ZERO;
        for (int k = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            suma = suma.add(get(k, i, j), mathContext);
        return suma.round(mathContext);
    }

    public final BigComplex promedio() {
        return sumaElementos().divide(BigComplex.valueOf(numeroElementos()), mathContext).round(mathContext);
    }

    public final P sumaPorCapas() {
        BigComplex[][] datos = new BigComplex[filas][columnas];
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++) {
            BigComplex suma = BigComplex.ZERO;
            for (int k = 0; k < capas; k++) suma = suma.add(get(k, i, j), mathContext);
            datos[i][j] = suma.round(mathContext);
        }
        return crearPlano(datos, mathContext, tolerancia);
    }

    public final P promedioPorCapas() { return sumaPorCapas().dividirPorEscalar(BigDecimal.valueOf(capas)); }

    public final BigComplex productoInternoFrobenius(Matriz3DCompleja otra) {
        mismasDimensiones(otra);
        BigComplex suma = BigComplex.ZERO;
        for (int k = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            suma = suma.add(get(k, i, j).conjugate().multiply(otra.get(k, i, j), mathContext), mathContext);
        return suma.round(mathContext);
    }

    public final BigComplex productoFrobenius(Matriz3DCompleja otra) {
        return productoInternoFrobenius(otra);
    }

    public final BigDecimal normaFrobenius() {
        return productoInternoFrobenius(this).re.max(BigDecimal.ZERO).sqrt(mathContext);
    }

    public final BigDecimal normaMaximo() {
        BigDecimal maximo = BigDecimal.ZERO;
        for (int k = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            maximo = maximo.max(get(k, i, j).abs(mathContext));
        return maximo;
    }

    /** Aplana en orden capa-fila-columna. */
    public final BigComplex[] flatten() {
        BigComplex[] salida = new BigComplex[numeroElementos()];
        for (int k = 0, n = 0; k < capas; k++)
            for (int i = 0; i < filas; i++)
                for (int j = 0; j < columnas; j++) salida[n++] = get(k, i, j);
        return salida;
    }

    public final T reshape(int nuevasCapas, int nuevasFilas, int nuevasColumnas) {
        int elementos;
        try { elementos = Math.multiplyExact(nuevasCapas, Math.multiplyExact(nuevasFilas, nuevasColumnas)); }
        catch (ArithmeticException e) { throw new DimensionMatrizException("Dimensiones demasiado grandes"); }
        if (nuevasCapas <= 0 || nuevasFilas <= 0 || nuevasColumnas <= 0 || elementos != numeroElementos())
            throw new DimensionMatrizException("reshape debe conservar " + numeroElementos() + " elementos");
        BigComplex[] plano = flatten();
        BigComplex[][][] datos = new BigComplex[nuevasCapas][nuevasFilas][nuevasColumnas];
        for (int k = 0, n = 0; k < nuevasCapas; k++)
            for (int i = 0; i < nuevasFilas; i++)
                for (int j = 0; j < nuevasColumnas; j++) datos[k][i][j] = plano[n++];
        return crear(datos, mathContext, tolerancia);
    }

    public final boolean esNula() {
        for (BigComplex valor : flatten()) if (!esCero(valor)) return false;
        return true;
    }

    public final boolean igualExacta(Matriz3DCompleja otra) {
        if (otra == null || capas != otra.capas() || filas != otra.filas() || columnas != otra.columnas()) return false;
        for (int k = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            if (!get(k, i, j).equals(otra.get(k, i, j))) return false;
        return true;
    }

    public final boolean aproximadamenteIgual(Matriz3DCompleja otra, BigDecimal tol) {
        Complejos.validarTolerancia(tol);
        if (otra == null || capas != otra.capas() || filas != otra.filas() || columnas != otra.columnas()) return false;
        for (int k = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            if (!Complejos.aproximadamenteIguales(get(k, i, j), otra.get(k, i, j), tol, mathContext)) return false;
        return true;
    }

    public final String formato(int escala) {
        if (escala < 0) throw new IllegalArgumentException("La escala no puede ser negativa");
        StringBuilder salida = new StringBuilder();
        for (int k = 0; k < capas; k++) {
            if (k > 0) salida.append(System.lineSeparator());
            salida.append("capa ").append(k).append(':').append(System.lineSeparator());
            for (int i = 0; i < filas; i++) {
                salida.append('[');
                for (int j = 0; j < columnas; j++) {
                    if (j > 0) salida.append(", ");
                    BigComplex z = get(k, i, j);
                    salida.append(z.re.setScale(escala, java.math.RoundingMode.HALF_UP).stripTrailingZeros());
                    if (z.im.signum() != 0) salida.append(z.im.signum() > 0 ? "+" : "")
                            .append(z.im.setScale(escala, java.math.RoundingMode.HALF_UP).stripTrailingZeros()).append('i');
                }
                salida.append(']').append(System.lineSeparator());
            }
        }
        return salida.toString().stripTrailing();
    }

    /** Serializa en orden capa-fila-columna sin pérdida decimal. */
    public final String serializar() {
        StringBuilder salida = new StringBuilder(dimensiones()).append(':');
        BigComplex[] valores = flatten();
        for (int n = 0; n < valores.length; n++) {
            if (n > 0) salida.append(';');
            salida.append(valores[n].re.toPlainString()).append(',').append(valores[n].im.toPlainString());
        }
        return salida.toString();
    }

    @Override public final boolean equals(Object objeto) {
        return objeto instanceof Matriz3DCompleja otra && igualExacta(otra);
    }

    @Override public final int hashCode() {
        int resultado = 1;
        for (BigComplex valor : flatten()) resultado = 31 * resultado + valor.hashCode();
        return 31 * (31 * (31 * resultado + capas) + filas) + columnas;
    }

    @Override public final String toString() { return formato(Math.min(8, mathContext.getPrecision())); }

    protected static BigComplex[][][] nuevosDatos(int capas, int filas, int columnas, BigComplex valor) {
        if (capas <= 0 || filas <= 0 || columnas <= 0) throw new DimensionMatrizException("Dimensiones 3D positivas requeridas");
        Objects.requireNonNull(valor, "El valor no puede ser null");
        BigComplex[][][] datos = new BigComplex[capas][filas][columnas];
        for (BigComplex[][] capa : datos) for (BigComplex[] fila : capa) Arrays.fill(fila, valor);
        return datos;
    }

    private P capaDe(Matriz3DCompleja origen, int capa) {
        BigComplex[][] datos = new BigComplex[origen.filas()][origen.columnas()];
        for (int i = 0; i < origen.filas(); i++)
            for (int j = 0; j < origen.columnas(); j++) datos[i][j] = origen.get(capa, i, j);
        return crearPlano(datos, mathContext, tolerancia);
    }

    private T ensamblarResultados(List<P> resultados) {
        P primera = resultados.getFirst();
        int nuevasFilas = primera.filas();
        int nuevasColumnas = primera.columnas();
        BigComplex[][][] datos = new BigComplex[capas][nuevasFilas][nuevasColumnas];
        for (int k = 0; k < capas; k++) {
            P resultado = resultados.get(k);
            if (!mathContext.equals(resultado.mathContext()) || tolerancia.compareTo(resultado.tolerancia()) != 0)
                throw new OperacionNoDefinidaException("Las capas resultantes cambiaron la política de precisión");
            if (resultado.filas() != nuevasFilas || resultado.columnas() != nuevasColumnas)
                throw new DimensionMatrizException("Las capas resultantes deben tener las mismas dimensiones");
            for (int i = 0; i < nuevasFilas; i++)
                for (int j = 0; j < nuevasColumnas; j++) datos[k][i][j] = resultado.get(i, j);
        }
        return crear(datos, mathContext, tolerancia);
    }

    private static void validarEje(int eje) {
        if (eje < 0 || eje > 2) throw new IllegalArgumentException("El eje debe ser 0, 1 o 2");
    }

    private static void validarPermutacion(int[] p) {
        if (Arrays.stream(p).anyMatch(x -> x < 0 || x > 2) || Arrays.stream(p).distinct().count() != 3)
            throw new IllegalArgumentException("Los ejes deben ser una permutación de 0, 1 y 2");
    }
}
