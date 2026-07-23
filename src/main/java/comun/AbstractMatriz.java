package comun;

import ch.obermuhlner.math.big.BigComplex;
import ch.obermuhlner.math.big.BigComplexMath;
import comun.excepciones.DimensionMatrizException;
import comun.excepciones.DivisionPorCeroException;
import comun.excepciones.IndiceMatrizException;
import comun.excepciones.MatrizNoCuadradaException;
import comun.excepciones.MatrizNoHermiticaException;
import comun.excepciones.MatrizSingularException;
import comun.excepciones.NoDefinidaPositivaException;
import comun.excepciones.OperacionNoDefinidaException;
import comun.resultados.ResultadoEliminacion;
import comun.resultados.ResultadoAutovalor;
import comun.resultados.ResultadoDiagonalizacion;
import comun.resultados.ResultadoIterativo;
import comun.resultados.ResultadoLU;
import comun.resultados.ResultadoQR;
import comun.resultados.ResultadoSVD;
import comun.resultados.ResultadoSistema;
import comun.resultados.TipoSolucion;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Objects;
import java.util.Random;
import java.util.Set;
import java.util.function.BinaryOperator;
import java.util.function.UnaryOperator;

import static comun.Complejos.ONE;
import static comun.Complejos.ZERO;

/**
 * Núcleo genérico inmutable de los algoritmos. Las subclases conservan su
 * almacenamiento exigido y únicamente exponen acceso por valor. Se rechazan
 * dimensiones nulas para evitar estados 0×n no representables por arrays;
 * todos los cálculos se redondean con el contexto de la matriz receptora y no
 * se permite combinar matrices con contextos distintos.
 *
 * @param <M> tipo concreto del resultado
 */
public abstract class AbstractMatriz<M extends AbstractMatriz<M>> implements MatrizCompleja {
    private final int filas;
    private final int columnas;
    private final MathContext mathContext;
    private final BigDecimal tolerancia;

    protected AbstractMatriz(int filas, int columnas, MathContext mathContext, BigDecimal tolerancia) {
        if (filas <= 0 || columnas <= 0) {
            throw new DimensionMatrizException("Las dimensiones deben ser positivas: " + filas + "x" + columnas);
        }
        this.filas = filas;
        this.columnas = columnas;
        this.mathContext = Objects.requireNonNull(mathContext, "El MathContext no puede ser null");
        Complejos.validarTolerancia(tolerancia);
        this.tolerancia = tolerancia;
    }

    protected abstract M crear(BigComplex[][] datos, MathContext contexto, BigDecimal tol);

    @Override public final int filas() { return filas; }
    @Override public final int columnas() { return columnas; }
    public final int getFilas() { return filas; }
    public final int getColumnas() { return columnas; }
    @Override public final MathContext mathContext() { return mathContext; }
    @Override public final BigDecimal tolerancia() { return tolerancia; }
    public final boolean esCuadrada() { return filas == columnas; }
    public final boolean esRectangular() { return filas != columnas; }
    public final boolean esMatrizFila() { return filas == 1; }
    public final boolean esMatrizColumna() { return columnas == 1; }
    public final BigComplex obtener(int fila, int columna) { return get(fila, columna); }

    protected final void validarIndice(int fila, int columna) {
        if (fila < 0 || fila >= filas || columna < 0 || columna >= columnas) {
            throw new IndiceMatrizException("Índice (" + fila + "," + columna + ") fuera de " + filas + "x" + columnas);
        }
    }

    protected final BigComplex r(BigComplex valor) {
        return valor.round(mathContext);
    }

    protected final boolean cero(BigComplex valor) {
        return Complejos.esCero(valor, tolerancia, mathContext);
    }

    protected final void mismaPolitica(MatrizCompleja otra) {
        Objects.requireNonNull(otra, "La matriz no puede ser null");
        if (!mathContext.equals(otra.mathContext()) || tolerancia.compareTo(otra.tolerancia()) != 0) {
            throw new OperacionNoDefinidaException("No se mezclan silenciosamente políticas de precisión distintas");
        }
    }

    protected final void mismasDimensiones(MatrizCompleja otra) {
        mismaPolitica(otra);
        if (filas != otra.filas() || columnas != otra.columnas()) {
            throw new DimensionMatrizException("Dimensiones incompatibles: " + filas + "x" + columnas
                    + " y " + otra.filas() + "x" + otra.columnas());
        }
    }

    protected final void exigirCuadrada() {
        if (!esCuadrada()) {
            throw new MatrizNoCuadradaException("La operación requiere matriz cuadrada, no " + filas + "x" + columnas);
        }
    }

    protected final BigComplex[][] copiaDatos() {
        BigComplex[][] copia = new BigComplex[filas][columnas];
        for (int i = 0; i < filas; i++) {
            for (int j = 0; j < columnas; j++) {
                copia[i][j] = get(i, j);
            }
        }
        return copia;
    }

    public final BigComplex[][] getMatriz() { return copiaDatos(); }

    public final ArrayList<ArrayList<BigComplex>> comoListaProfunda() {
        ArrayList<ArrayList<BigComplex>> salida = new ArrayList<>(filas);
        for (int i = 0; i < filas; i++) {
            ArrayList<BigComplex> fila = new ArrayList<>(columnas);
            Collections.addAll(fila, obtenerFila(i));
            salida.add(fila);
        }
        return salida;
    }

    public final array.Matriz comoMatrizArray() {
        return new array.Matriz(this, mathContext, tolerancia);
    }

    public final lista.Matriz comoMatrizLista() {
        return new lista.Matriz(this, mathContext, tolerancia);
    }

    public final BigComplex[] obtenerFila(int i) {
        validarIndice(i, 0);
        BigComplex[] salida = new BigComplex[columnas];
        for (int j = 0; j < columnas; j++) salida[j] = get(i, j);
        return salida;
    }

    public final BigComplex[] obtenerColumna(int j) {
        validarIndice(0, j);
        BigComplex[] salida = new BigComplex[filas];
        for (int i = 0; i < filas; i++) salida[i] = get(i, j);
        return salida;
    }

    public final M conElemento(int fila, int columna, BigComplex valor) {
        validarIndice(fila, columna);
        BigComplex[][] d = copiaDatos();
        d[fila][columna] = r(Objects.requireNonNull(valor, "El valor no puede ser null"));
        return crear(d, mathContext, tolerancia);
    }

    public final M conFila(int fila, BigComplex[] valores) {
        validarIndice(fila, 0);
        Objects.requireNonNull(valores, "La fila no puede ser null");
        if (valores.length != columnas) throw new DimensionMatrizException("Longitud de fila incorrecta");
        BigComplex[][] d = copiaDatos();
        for (int j = 0; j < columnas; j++) d[fila][j] = r(Objects.requireNonNull(valores[j]));
        return crear(d, mathContext, tolerancia);
    }

    public final M conColumna(int columna, BigComplex[] valores) {
        validarIndice(0, columna);
        Objects.requireNonNull(valores, "La columna no puede ser null");
        if (valores.length != filas) throw new DimensionMatrizException("Longitud de columna incorrecta");
        BigComplex[][] d = copiaDatos();
        for (int i = 0; i < filas; i++) d[i][columna] = r(Objects.requireNonNull(valores[i]));
        return crear(d, mathContext, tolerancia);
    }

    /** Submatriz con límites desde inclusivos y hasta exclusivos. */
    public final M submatriz(int filaDesde, int filaHasta, int columnaDesde, int columnaHasta) {
        if (filaDesde < 0 || filaHasta > filas || columnaDesde < 0 || columnaHasta > columnas
                || filaDesde >= filaHasta || columnaDesde >= columnaHasta) {
            throw new IndiceMatrizException("Intervalo de submatriz inválido");
        }
        BigComplex[][] d = new BigComplex[filaHasta - filaDesde][columnaHasta - columnaDesde];
        for (int i = filaDesde; i < filaHasta; i++)
            for (int j = columnaDesde; j < columnaHasta; j++)
                d[i - filaDesde][j - columnaDesde] = get(i, j);
        return crear(d, mathContext, tolerancia);
    }

    public final M extraerFilas(int... indices) {
        Objects.requireNonNull(indices);
        if (indices.length == 0) throw new DimensionMatrizException("Debe extraerse alguna fila");
        BigComplex[][] d = new BigComplex[indices.length][columnas];
        for (int i = 0; i < indices.length; i++) {
            validarIndice(indices[i], 0);
            d[i] = obtenerFila(indices[i]);
        }
        return crear(d, mathContext, tolerancia);
    }

    public final M extraerColumnas(int... indices) {
        Objects.requireNonNull(indices);
        if (indices.length == 0) throw new DimensionMatrizException("Debe extraerse alguna columna");
        BigComplex[][] d = new BigComplex[filas][indices.length];
        for (int j = 0; j < indices.length; j++) {
            validarIndice(0, indices[j]);
            for (int i = 0; i < filas; i++) d[i][j] = get(i, indices[j]);
        }
        return crear(d, mathContext, tolerancia);
    }

    public final M eliminarFila(int fila) {
        if (filas == 1) throw new DimensionMatrizException("No se admite una matriz sin filas");
        validarIndice(fila, 0);
        int[] idx = new int[filas - 1];
        for (int i = 0, k = 0; i < filas; i++) if (i != fila) idx[k++] = i;
        return extraerFilas(idx);
    }

    public final M eliminarColumna(int columna) {
        if (columnas == 1) throw new DimensionMatrizException("No se admite una matriz sin columnas");
        validarIndice(0, columna);
        int[] idx = new int[columnas - 1];
        for (int j = 0, k = 0; j < columnas; j++) if (j != columna) idx[k++] = j;
        return extraerColumnas(idx);
    }

    public final M eliminarFilaYColumna(int fila, int columna) {
        return eliminarFila(fila).eliminarColumna(columna);
    }

    public final M concatenarHorizontal(MatrizCompleja otra) {
        mismaPolitica(otra);
        if (filas != otra.filas()) throw new DimensionMatrizException("Las filas deben coincidir");
        BigComplex[][] d = new BigComplex[filas][columnas + otra.columnas()];
        for (int i = 0; i < filas; i++) {
            for (int j = 0; j < columnas; j++) d[i][j] = get(i, j);
            for (int j = 0; j < otra.columnas(); j++) d[i][columnas + j] = otra.get(i, j);
        }
        return crear(d, mathContext, tolerancia);
    }

    public final M concatenarVertical(MatrizCompleja otra) {
        mismaPolitica(otra);
        if (columnas != otra.columnas()) throw new DimensionMatrizException("Las columnas deben coincidir");
        BigComplex[][] d = new BigComplex[filas + otra.filas()][columnas];
        for (int i = 0; i < filas; i++) d[i] = obtenerFila(i);
        for (int i = 0; i < otra.filas(); i++)
            for (int j = 0; j < columnas; j++) d[filas + i][j] = otra.get(i, j);
        return crear(d, mathContext, tolerancia);
    }

    public final M insertarFila(int posicion, BigComplex[] valores) {
        if (posicion < 0 || posicion > filas) throw new IndiceMatrizException("Posición de fila inválida");
        if (valores == null || valores.length != columnas) throw new DimensionMatrizException("Fila incompatible");
        BigComplex[][] d = new BigComplex[filas + 1][columnas];
        for (int i = 0, k = 0; i < d.length; i++) {
            if (i == posicion) {
                for (int j = 0; j < columnas; j++) d[i][j] = r(Objects.requireNonNull(valores[j]));
            } else d[i] = obtenerFila(k++);
        }
        return crear(d, mathContext, tolerancia);
    }

    public final M insertarColumna(int posicion, BigComplex[] valores) {
        if (posicion < 0 || posicion > columnas) throw new IndiceMatrizException("Posición de columna inválida");
        if (valores == null || valores.length != filas) throw new DimensionMatrizException("Columna incompatible");
        BigComplex[][] d = new BigComplex[filas][columnas + 1];
        for (int i = 0; i < filas; i++) {
            for (int j = 0, k = 0; j < d[i].length; j++) {
                d[i][j] = j == posicion ? r(Objects.requireNonNull(valores[i])) : get(i, k++);
            }
        }
        return crear(d, mathContext, tolerancia);
    }

    public final M redimensionar(int nuevasFilas, int nuevasColumnas, BigComplex relleno) {
        if (nuevasFilas <= 0 || nuevasColumnas <= 0) throw new DimensionMatrizException("Dimensión inválida");
        BigComplex v = r(Objects.requireNonNull(relleno));
        BigComplex[][] d = nuevosDatos(nuevasFilas, nuevasColumnas, v);
        for (int i = 0; i < Math.min(filas, nuevasFilas); i++)
            for (int j = 0; j < Math.min(columnas, nuevasColumnas); j++) d[i][j] = get(i, j);
        return crear(d, mathContext, tolerancia);
    }

    public final M rellenar(BigComplex valor) {
        return crear(nuevosDatos(filas, columnas, r(Objects.requireNonNull(valor))), mathContext, tolerancia);
    }

    public final M map(UnaryOperator<BigComplex> funcion) {
        Objects.requireNonNull(funcion);
        BigComplex[][] d = new BigComplex[filas][columnas];
        for (int i = 0; i < filas; i++)
            for (int j = 0; j < columnas; j++) d[i][j] = r(Objects.requireNonNull(funcion.apply(get(i, j))));
        return crear(d, mathContext, tolerancia);
    }

    public final M zip(MatrizCompleja otra, BinaryOperator<BigComplex> funcion) {
        mismasDimensiones(otra);
        Objects.requireNonNull(funcion);
        BigComplex[][] d = new BigComplex[filas][columnas];
        for (int i = 0; i < filas; i++)
            for (int j = 0; j < columnas; j++)
                d[i][j] = r(Objects.requireNonNull(funcion.apply(get(i, j), otra.get(i, j))));
        return crear(d, mathContext, tolerancia);
    }

    private static BigComplex[][] nuevosDatos(int f, int c, BigComplex valor) {
        BigComplex[][] d = new BigComplex[f][c];
        for (BigComplex[] fila : d) Arrays.fill(fila, valor);
        return d;
    }

    public final M suma(MatrizCompleja otra) {
        return zip(otra, (a, b) -> a.add(b, mathContext));
    }

    public final M resta(MatrizCompleja otra) {
        return zip(otra, (a, b) -> a.subtract(b, mathContext));
    }

    public final M opuesta() { return map(BigComplex::negate); }

    public final M productoEscalar(BigComplex escalar) {
        Objects.requireNonNull(escalar);
        return map(x -> x.multiply(escalar, mathContext));
    }

    public final M productoEscalar(BigDecimal escalar) {
        return productoEscalar(BigComplex.valueOf(Objects.requireNonNull(escalar)));
    }

    public final M productoEscalar(double escalar) {
        return productoEscalar(Complejos.de(escalar));
    }

    public final M dividirPorEscalar(BigComplex escalar) {
        Objects.requireNonNull(escalar);
        if (cero(escalar)) throw new DivisionPorCeroException("No se puede dividir por cero");
        return map(x -> x.divide(escalar, mathContext));
    }

    public final M dividirPorEscalar(BigDecimal escalar) {
        return dividirPorEscalar(BigComplex.valueOf(Objects.requireNonNull(escalar)));
    }

    public final M combinacionLineal(BigComplex alfa, MatrizCompleja otra, BigComplex beta) {
        return productoEscalar(alfa).suma(comoEste(otra).productoEscalar(beta));
    }

    @SuppressWarnings("unchecked")
    private M comoEste(MatrizCompleja otra) {
        mismaPolitica(otra);
        if (getClass().isInstance(otra)) return (M) otra;
        BigComplex[][] d = new BigComplex[otra.filas()][otra.columnas()];
        for (int i = 0; i < otra.filas(); i++)
            for (int j = 0; j < otra.columnas(); j++) d[i][j] = otra.get(i, j);
        return crear(d, mathContext, tolerancia);
    }

    public final M producto(MatrizCompleja otra) {
        mismaPolitica(otra);
        if (columnas != otra.filas()) throw new DimensionMatrizException("Producto matricial incompatible");
        BigComplex[][] d = nuevosDatos(filas, otra.columnas(), ZERO);
        for (int i = 0; i < filas; i++)
            for (int k = 0; k < columnas; k++)
                for (int j = 0; j < otra.columnas(); j++)
                    d[i][j] = r(d[i][j].add(get(i, k).multiply(otra.get(k, j), mathContext), mathContext));
        return crear(d, mathContext, tolerancia);
    }

    public final BigComplex[] productoPorVectorColumna(BigComplex[] vector) {
        Objects.requireNonNull(vector);
        if (vector.length != columnas) throw new DimensionMatrizException("Vector incompatible");
        BigComplex[] salida = new BigComplex[filas];
        for (int i = 0; i < filas; i++) {
            BigComplex suma = ZERO;
            for (int j = 0; j < columnas; j++)
                suma = suma.add(get(i, j).multiply(Objects.requireNonNull(vector[j]), mathContext), mathContext);
            salida[i] = r(suma);
        }
        return salida;
    }

    public final BigComplex[] vectorFilaPorMatriz(BigComplex[] vector) {
        Objects.requireNonNull(vector);
        if (vector.length != filas) throw new DimensionMatrizException("Vector incompatible");
        BigComplex[] salida = new BigComplex[columnas];
        for (int j = 0; j < columnas; j++) {
            BigComplex suma = ZERO;
            for (int i = 0; i < filas; i++)
                suma = suma.add(Objects.requireNonNull(vector[i]).multiply(get(i, j), mathContext), mathContext);
            salida[j] = r(suma);
        }
        return salida;
    }

    public final M potencia(int exponente) { return potenciaRapida(exponente); }

    public final M potenciaRapida(int exponente) {
        exigirCuadrada();
        if (exponente < 0) return potenciaNegativa(exponente);
        M resultado = identidadInterna(filas);
        M base = comoEste(this);
        int e = exponente;
        while (e > 0) {
            if ((e & 1) == 1) resultado = resultado.producto(base);
            e >>>= 1;
            if (e > 0) base = base.producto(base);
        }
        return resultado;
    }

    public final M potenciaNegativa(int exponente) {
        if (exponente >= 0) throw new IllegalArgumentException("El exponente debe ser negativo");
        return inversa().potenciaRapida(-exponente);
    }

    public final M productoHadamard(MatrizCompleja otra) {
        return zip(otra, (a, b) -> a.multiply(b, mathContext));
    }

    public final M cocienteHadamard(MatrizCompleja otra) {
        return zip(otra, (a, b) -> {
            if (cero(b)) throw new DivisionPorCeroException("Cero en el divisor de Hadamard");
            return a.divide(b, mathContext);
        });
    }

    public final M productoKronecker(MatrizCompleja otra) {
        mismaPolitica(otra);
        BigComplex[][] d = new BigComplex[filas * otra.filas()][columnas * otra.columnas()];
        for (int i = 0; i < filas; i++)
            for (int j = 0; j < columnas; j++)
                for (int k = 0; k < otra.filas(); k++)
                    for (int l = 0; l < otra.columnas(); l++)
                        d[i * otra.filas() + k][j * otra.columnas() + l] =
                                r(get(i, j).multiply(otra.get(k, l), mathContext));
        return crear(d, mathContext, tolerancia);
    }

    public final M productoKhatriRao(MatrizCompleja otra) {
        mismaPolitica(otra);
        if (columnas != otra.columnas()) throw new DimensionMatrizException("Khatri-Rao exige igual número de columnas");
        BigComplex[][] d = new BigComplex[filas * otra.filas()][columnas];
        for (int j = 0; j < columnas; j++)
            for (int i = 0; i < filas; i++)
                for (int k = 0; k < otra.filas(); k++)
                    d[i * otra.filas() + k][j] = r(get(i, j).multiply(otra.get(k, j), mathContext));
        return crear(d, mathContext, tolerancia);
    }

    public final M productoExterior(BigComplex[] x, BigComplex[] y) {
        Objects.requireNonNull(x); Objects.requireNonNull(y);
        BigComplex[][] d = new BigComplex[x.length][y.length];
        for (int i = 0; i < x.length; i++)
            for (int j = 0; j < y.length; j++)
                d[i][j] = r(Objects.requireNonNull(x[i]).multiply(Objects.requireNonNull(y[j]).conjugate(), mathContext));
        return crear(d, mathContext, tolerancia);
    }

    public final BigComplex productoInterior(BigComplex[] x, BigComplex[] y) {
        if (x == null || y == null || x.length != y.length) throw new DimensionMatrizException("Vectores incompatibles");
        BigComplex s = ZERO;
        for (int i = 0; i < x.length; i++)
            s = s.add(Objects.requireNonNull(x[i]).conjugate().multiply(Objects.requireNonNull(y[i]), mathContext), mathContext);
        return r(s);
    }

    public final BigComplex productoFrobenius(MatrizCompleja otra) {
        mismasDimensiones(otra);
        BigComplex s = ZERO;
        for (int i = 0; i < filas; i++)
            for (int j = 0; j < columnas; j++)
                s = s.add(get(i, j).conjugate().multiply(otra.get(i, j), mathContext), mathContext);
        return r(s);
    }

    public final BigComplex productoInternoFrobenius(MatrizCompleja otra) { return productoFrobenius(otra); }
    public final M conmutador(MatrizCompleja otra) { return producto(otra).resta(comoEste(otra).producto(this)); }
    public final M anticonmutador(MatrizCompleja otra) { return producto(otra).suma(comoEste(otra).producto(this)); }

    public final BigComplex formaBilineal(BigComplex[] x, BigComplex[] y) {
        return productoInteriorSinConjugar(x, productoPorVectorColumna(y));
    }

    public final BigComplex formaSesquilineal(BigComplex[] x, BigComplex[] y) {
        return productoInterior(x, productoPorVectorColumna(y));
    }

    public final BigComplex formaCuadratica(BigComplex[] x) { return formaSesquilineal(x, x); }
    public final BigComplex quadraticForm(BigComplex[] x) { return formaCuadratica(x); }

    public final BigComplex sandwich(BigComplex[] x, BigComplex[] y) {
        return formaSesquilineal(x, y);
    }

    private BigComplex productoInteriorSinConjugar(BigComplex[] x, BigComplex[] y) {
        if (x == null || y == null || x.length != y.length) throw new DimensionMatrizException("Vectores incompatibles");
        BigComplex s = ZERO;
        for (int i = 0; i < x.length; i++) s = s.add(x[i].multiply(y[i], mathContext), mathContext);
        return r(s);
    }

    public final M traspuesta() {
        BigComplex[][] d = new BigComplex[columnas][filas];
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++) d[j][i] = get(i, j);
        return crear(d, mathContext, tolerancia);
    }

    public final M conjugada() { return map(BigComplex::conjugate); }
    public final M adjuntaHermitica() { return conjugada().traspuesta(); }
    public final M transpuestaConjugada() { return adjuntaHermitica(); }
    public final M parteReal() { return map(z -> BigComplex.valueOf(z.re)); }
    public final M parteImaginaria() { return map(z -> BigComplex.valueOf(z.im)); }
    public final M parteSimetrica() { exigirCuadrada(); return suma(traspuesta()).dividirPorEscalar(BigDecimal.TWO); }
    public final M parteAntisimetrica() { exigirCuadrada(); return resta(traspuesta()).dividirPorEscalar(BigDecimal.TWO); }
    public final M parteHermitica() { exigirCuadrada(); return suma(adjuntaHermitica()).dividirPorEscalar(BigDecimal.TWO); }
    public final M parteAntihermitica() { exigirCuadrada(); return resta(adjuntaHermitica()).dividirPorEscalar(BigDecimal.TWO); }

    public final BigComplex[] diagonalPrincipal() { return diagonales(0); }

    public final BigComplex[] diagonales(int desplazamiento) {
        if (desplazamiento <= -filas || desplazamiento >= columnas)
            throw new IndiceMatrizException("Diagonal inexistente");
        int i0 = Math.max(0, -desplazamiento), j0 = Math.max(0, desplazamiento);
        int n = Math.min(filas - i0, columnas - j0);
        BigComplex[] d = new BigComplex[n];
        for (int k = 0; k < n; k++) d[k] = get(i0 + k, j0 + k);
        return d;
    }

    public final M triangularSuperior() {
        return mapPosicional((i, j, z) -> i <= j ? z : ZERO);
    }

    public final M triangularInferior() {
        return mapPosicional((i, j, z) -> i >= j ? z : ZERO);
    }

    public final M estrictamenteTriangularSuperior() {
        return mapPosicional((i, j, z) -> i < j ? z : ZERO);
    }

    public final M estrictamenteTriangularInferior() {
        return mapPosicional((i, j, z) -> i > j ? z : ZERO);
    }

    @FunctionalInterface
    protected interface OperacionPosicional {
        BigComplex aplicar(int fila, int columna, BigComplex valor);
    }

    protected final M mapPosicional(OperacionPosicional op) {
        BigComplex[][] d = new BigComplex[filas][columnas];
        for (int i = 0; i < filas; i++)
            for (int j = 0; j < columnas; j++) d[i][j] = op.aplicar(i, j, get(i, j));
        return crear(d, mathContext, tolerancia);
    }

    protected final M identidadInterna(int n) {
        BigComplex[][] d = nuevosDatos(n, n, ZERO);
        for (int i = 0; i < n; i++) d[i][i] = ONE;
        return crear(d, mathContext, tolerancia);
    }

    public final boolean esNula() {
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++) if (!cero(get(i, j))) return false;
        return true;
    }

    public final boolean esIdentidad() { return esCuadrada() && esIdentidadRectangular(); }

    public final boolean esIdentidadRectangular() {
        for (int i = 0; i < filas; i++)
            for (int j = 0; j < columnas; j++) {
                BigComplex esperado = i == j ? ONE : ZERO;
                if (!Complejos.aproximadamenteIguales(get(i, j), esperado, tolerancia, mathContext)) return false;
            }
        return true;
    }

    public final boolean esDiagonal() {
        for (int i = 0; i < filas; i++)
            for (int j = 0; j < columnas; j++) if (i != j && !cero(get(i, j))) return false;
        return true;
    }

    public final boolean esEscalar() {
        if (!esCuadrada() || !esDiagonal()) return false;
        for (int i = 1; i < filas; i++)
            if (!Complejos.aproximadamenteIguales(get(0, 0), get(i, i), tolerancia, mathContext)) return false;
        return true;
    }

    public final boolean esTriangularSuperior() {
        for (int i = 1; i < filas; i++) for (int j = 0; j < Math.min(i, columnas); j++) if (!cero(get(i, j))) return false;
        return true;
    }

    public final boolean esTriangularInferior() {
        for (int i = 0; i < filas; i++) for (int j = i + 1; j < columnas; j++) if (!cero(get(i, j))) return false;
        return true;
    }

    public final boolean esEstrictamenteTriangularSuperior() {
        return esTriangularSuperior() && diagonalCero();
    }

    public final boolean esEstrictamenteTriangularInferior() {
        return esTriangularInferior() && diagonalCero();
    }

    private boolean diagonalCero() {
        for (int i = 0; i < Math.min(filas, columnas); i++) if (!cero(get(i, i))) return false;
        return true;
    }

    public final boolean esBidiagonalSuperior() {
        return soloBandas(Set.of(0, 1));
    }

    public final boolean esBidiagonalInferior() {
        return soloBandas(Set.of(0, -1));
    }

    public final boolean esTridiagonal() {
        return soloBandas(Set.of(-1, 0, 1));
    }

    private boolean soloBandas(Set<Integer> bandas) {
        for (int i = 0; i < filas; i++)
            for (int j = 0; j < columnas; j++) if (!bandas.contains(j - i) && !cero(get(i, j))) return false;
        return true;
    }

    public final boolean esMatrizBanda() {
        return anchoDeBandaSuperior() + anchoDeBandaInferior() + 1 < Math.max(filas, columnas);
    }

    public final int anchoDeBandaSuperior() {
        int ancho = 0;
        for (int i = 0; i < filas; i++) for (int j = i; j < columnas; j++) if (!cero(get(i, j))) ancho = Math.max(ancho, j - i);
        return ancho;
    }

    public final int anchoDeBandaInferior() {
        int ancho = 0;
        for (int i = 0; i < filas; i++) for (int j = 0; j <= Math.min(i, columnas - 1); j++)
            if (!cero(get(i, j))) ancho = Math.max(ancho, i - j);
        return ancho;
    }

    public final boolean esHessenbergSuperior() { return anchoDeBandaInferior() <= 1; }
    public final boolean esHessenbergInferior() { return anchoDeBandaSuperior() <= 1; }

    public final boolean esToeplitz() {
        for (int i = 1; i < filas; i++) for (int j = 1; j < columnas; j++)
            if (!Complejos.aproximadamenteIguales(get(i, j), get(i - 1, j - 1), tolerancia, mathContext)) return false;
        return true;
    }

    public final boolean esHankel() {
        for (int i = 1; i < filas; i++) for (int j = 0; j + 1 < columnas; j++)
            if (!Complejos.aproximadamenteIguales(get(i, j), get(i - 1, j + 1), tolerancia, mathContext)) return false;
        return true;
    }

    public final boolean esCirculante() {
        if (!esCuadrada()) return false;
        for (int i = 1; i < filas; i++) for (int j = 0; j < columnas; j++)
            if (!Complejos.aproximadamenteIguales(get(i, j), get(i - 1, (j - 1 + columnas) % columnas),
                    tolerancia, mathContext)) return false;
        return true;
    }

    public final boolean esMatrizDePermutacion() {
        if (!esCuadrada()) return false;
        for (int i = 0; i < filas; i++) {
            int unos = 0;
            for (int j = 0; j < columnas; j++) {
                if (Complejos.aproximadamenteIguales(get(i, j), ONE, tolerancia, mathContext)) unos++;
                else if (!cero(get(i, j))) return false;
            }
            if (unos != 1) return false;
        }
        for (int j = 0; j < columnas; j++) {
            int unos = 0;
            for (int i = 0; i < filas; i++) if (Complejos.aproximadamenteIguales(get(i, j), ONE, tolerancia, mathContext)) unos++;
            if (unos != 1) return false;
        }
        return true;
    }

    public final boolean esMonomial() {
        if (!esCuadrada()) return false;
        for (int i = 0; i < filas; i++) {
            int nz = 0;
            for (int j = 0; j < columnas; j++) if (!cero(get(i, j))) nz++;
            if (nz != 1) return false;
        }
        for (int j = 0; j < columnas; j++) {
            int nz = 0;
            for (int i = 0; i < filas; i++) if (!cero(get(i, j))) nz++;
            if (nz != 1) return false;
        }
        return true;
    }

    public final boolean esMatrizPorBloques(int[] cortesFilas, int[] cortesColumnas) {
        validarParticion(cortesFilas, filas);
        validarParticion(cortesColumnas, columnas);
        return true;
    }

    public final boolean esDiagonalPorBloques(int[] tamanos) {
        exigirCuadrada();
        validarTamanosBloques(tamanos, filas);
        int inicioI = 0;
        for (int bi = 0; bi < tamanos.length; bi++) {
            int inicioJ = 0;
            for (int bj = 0; bj < tamanos.length; bj++) {
                if (bi != bj) {
                    for (int i = inicioI; i < inicioI + tamanos[bi]; i++)
                        for (int j = inicioJ; j < inicioJ + tamanos[bj]; j++) if (!cero(get(i, j))) return false;
                }
                inicioJ += tamanos[bj];
            }
            inicioI += tamanos[bi];
        }
        return true;
    }

    public final BigDecimal densidad() {
        int nz = 0;
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++) if (!cero(get(i, j))) nz++;
        return BigDecimal.valueOf(nz).divide(BigDecimal.valueOf((long) filas * columnas), mathContext);
    }

    public final boolean esDispersa() { return densidad().compareTo(new BigDecimal("0.2")) <= 0; }

    public final boolean esBinaria() {
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            if (!cero(get(i, j)) && !Complejos.aproximadamenteIguales(get(i, j), ONE, tolerancia, mathContext)) return false;
        return true;
    }

    public final boolean esMatrizDeSignos() {
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++) {
            BigComplex z = get(i, j);
            if (!z.isReal() || (z.re.compareTo(BigDecimal.ONE) != 0 && z.re.compareTo(BigDecimal.ONE.negate()) != 0)) return false;
        }
        return true;
    }

    public final boolean esNoNegativa() {
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            if (!get(i, j).isReal() || get(i, j).re.compareTo(tolerancia.negate()) < 0) return false;
        return true;
    }

    public final boolean esPositivaPorEntradas() {
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            if (!get(i, j).isReal() || get(i, j).re.compareTo(tolerancia) <= 0) return false;
        return true;
    }

    public final boolean esEstocasticaPorFilas() {
        if (!esNoNegativa()) return false;
        for (int i = 0; i < filas; i++) {
            BigComplex s = ZERO;
            for (int j = 0; j < columnas; j++) s = s.add(get(i, j), mathContext);
            if (!Complejos.aproximadamenteIguales(s, ONE, tolerancia, mathContext)) return false;
        }
        return true;
    }

    public final boolean esEstocasticaPorColumnas() {
        return traspuesta().esEstocasticaPorFilas();
    }

    public final boolean esDoblementeEstocastica() { return esEstocasticaPorFilas() && esEstocasticaPorColumnas(); }
    public final boolean esMatrizDeTransicion() { return esCuadrada() && esEstocasticaPorFilas(); }
    public final boolean esMatrizTransicion() { return esMatrizDeTransicion(); }

    public final boolean esSimetrica() { return esCuadrada() && aproximadamenteIgual(traspuesta(), tolerancia); }
    public final boolean esAntisimetrica() { return esCuadrada() && aproximadamenteIgual(traspuesta().opuesta(), tolerancia); }
    public final boolean esHermitica() { return esCuadrada() && aproximadamenteIgual(adjuntaHermitica(), tolerancia); }
    public final boolean esAntihermitica() { return esCuadrada() && aproximadamenteIgual(adjuntaHermitica().opuesta(), tolerancia); }
    public final boolean esOrtogonal() { return esCuadrada() && traspuesta().producto(this).esIdentidad(); }
    public final boolean esUnitaria() { return esCuadrada() && adjuntaHermitica().producto(this).esIdentidad(); }
    public final boolean esNormal() { return esCuadrada() && conmutaCon(adjuntaHermitica()); }
    public final boolean esIdempotente() { return esCuadrada() && potencia(2).aproximadamenteIgual(this, tolerancia); }
    public final boolean esInvolutiva() { return esCuadrada() && potencia(2).esIdentidad(); }

    public final boolean esNilpotente() {
        if (!esCuadrada()) return false;
        M p = comoEste(this);
        for (int k = 1; k <= filas; k++) {
            if (p.esNula()) return true;
            p = p.producto(this);
        }
        return false;
    }

    public final int indiceDeNilpotencia() {
        if (!esCuadrada()) throw new MatrizNoCuadradaException("Nilpotencia exige matriz cuadrada");
        M p = comoEste(this);
        for (int k = 1; k <= filas; k++) {
            if (p.esNula()) return k;
            p = p.producto(this);
        }
        throw new OperacionNoDefinidaException("La matriz no es nilpotente");
    }

    public final boolean esPeriodica() {
        if (!esCuadrada()) return false;
        M p = comoEste(this);
        for (int k = 1; k <= 2 * filas * filas + 1; k++) {
            if (k > 1 && p.aproximadamenteIgual(this, tolerancia)) return true;
            p = p.producto(this);
        }
        return false;
    }

    public final int periodo() {
        exigirCuadrada();
        M p = comoEste(this);
        for (int k = 1; k <= 2 * filas * filas + 1; k++) {
            if (k > 1 && p.aproximadamenteIgual(this, tolerancia)) return k - 1;
            p = p.producto(this);
        }
        throw new OperacionNoDefinidaException("No se detectó periodicidad en la cota finita");
    }

    public final boolean esSingular() { return esCuadrada() && cero(determinante()); }
    public final boolean esRegular() { return esCuadrada() && !esSingular(); }
    public final boolean esInvertible() { return esRegular(); }
    public final boolean tieneRangoCompleto() { return rango() == Math.min(filas, columnas); }
    public final boolean tieneRangoCompletoPorFilas() { return rango() == filas; }
    public final boolean tieneRangoCompletoPorColumnas() { return rango() == columnas; }
    public final boolean esProyector() { return esIdempotente(); }
    public final boolean esProyectorOrtogonal() { return esProyector() && esHermitica(); }
    public final boolean esProyectorOblicuo() { return esProyector() && !esHermitica(); }

    public final boolean esRotacion() {
        return esCuadrada() && todosReales() && esOrtogonal()
                && Complejos.aproximadamenteIguales(determinante(), ONE, tolerancia, mathContext);
    }

    public final boolean esReflexion() {
        return esCuadrada() && todosReales() && esOrtogonal()
                && Complejos.aproximadamenteIguales(determinante(), Complejos.MINUS_ONE, tolerancia, mathContext);
    }

    public final boolean esReducible() {
        if (!esCuadrada()) return false;
        return componentesConexas().size() > 1;
    }

    public final boolean esIrreducible() { return esCuadrada() && !esReducible(); }
    public final boolean conmutaCon(MatrizCompleja otra) {
        if (!esCuadrada() || filas != otra.filas() || otra.filas() != otra.columnas()) return false;
        return producto(otra).aproximadamenteIgual(comoEste(otra).producto(this), tolerancia);
    }

    public final boolean esSemejanteA(MatrizCompleja otra) {
        if (!esCuadrada() || filas != otra.filas() || otra.filas() != otra.columnas()) return false;
        return Complejos.aproximadamenteIguales(traza(), comoEste(otra).traza(), tolerancia, mathContext)
                && Complejos.aproximadamenteIguales(determinante(), comoEste(otra).determinante(), tolerancia, mathContext)
                && rango() == comoEste(otra).rango();
    }

    /** Para hermíticas, la congruencia se decide por la inercia de Sylvester. */
    public final boolean esCongruenteCon(MatrizCompleja otra) {
        M b = comoEste(otra);
        if (!esHermitica() || !b.esHermitica() || filas != b.filas()) return false;
        return Arrays.equals(inercia(), b.inercia());
    }

    public final boolean esDiagonalmenteDominantePorFilas() {
        if (!esCuadrada()) return false;
        for (int i = 0; i < filas; i++) {
            BigDecimal suma = BigDecimal.ZERO;
            for (int j = 0; j < columnas; j++) if (i != j) suma = suma.add(get(i, j).abs(mathContext), mathContext);
            if (get(i, i).abs(mathContext).compareTo(suma) < 0) return false;
        }
        return true;
    }

    public final boolean esDiagonalmenteDominantePorColumnas() { return traspuesta().esDiagonalmenteDominantePorFilas(); }

    public final boolean esEstrictamenteDiagonalmenteDominante() {
        if (!esCuadrada()) return false;
        for (int i = 0; i < filas; i++) {
            BigDecimal suma = BigDecimal.ZERO;
            for (int j = 0; j < columnas; j++) if (i != j) suma = suma.add(get(i, j).abs(mathContext), mathContext);
            if (get(i, i).abs(mathContext).compareTo(suma) <= 0) return false;
        }
        return true;
    }

    public final boolean esHermiticaDefinidaPositiva() { return esHermitica() && esDefinidaPositiva(); }
    public final boolean esNormalAproximadamente(BigDecimal tol) {
        Complejos.validarTolerancia(tol);
        return esCuadrada() && producto(adjuntaHermitica()).aproximadamenteIgual(adjuntaHermitica().producto(this), tol);
    }

    public final boolean esUnitariaAproximadamente(BigDecimal tol) {
        Complejos.validarTolerancia(tol);
        return esCuadrada() && adjuntaHermitica().producto(this).aproximadamenteIgual(identidadInterna(filas), tol);
    }

    public final boolean igualExacta(MatrizCompleja otra) {
        if (otra == null || filas != otra.filas() || columnas != otra.columnas()) return false;
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            if (!get(i, j).equals(otra.get(i, j))) return false;
        return true;
    }

    public final boolean aproximadamenteIgual(MatrizCompleja otra, BigDecimal tol) {
        Complejos.validarTolerancia(tol);
        if (otra == null || filas != otra.filas() || columnas != otra.columnas()) return false;
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            if (!Complejos.aproximadamenteIguales(get(i, j), otra.get(i, j), tol, mathContext)) return false;
        return true;
    }

    private boolean todosReales() {
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++) if (!cero(BigComplex.valueOf(get(i, j).im))) return false;
        return true;
    }

    public final M intercambiarFilas(int i, int j) {
        validarIndice(i, 0); validarIndice(j, 0);
        BigComplex[][] d = copiaDatos();
        BigComplex[] tmp = d[i]; d[i] = d[j]; d[j] = tmp;
        return crear(d, mathContext, tolerancia);
    }

    public final M multiplicarFila(int i, BigComplex escalarNoNulo) {
        validarIndice(i, 0);
        if (escalarNoNulo == null || cero(escalarNoNulo))
            throw new DivisionPorCeroException("La operación elemental debe ser reversible");
        BigComplex[][] d = copiaDatos();
        for (int j = 0; j < columnas; j++) d[i][j] = r(d[i][j].multiply(escalarNoNulo, mathContext));
        return crear(d, mathContext, tolerancia);
    }

    public final M sumarMultiploDeFila(int origen, int destino, BigComplex escalar) {
        validarIndice(origen, 0); validarIndice(destino, 0); Objects.requireNonNull(escalar);
        BigComplex[][] d = copiaDatos();
        for (int j = 0; j < columnas; j++)
            d[destino][j] = r(d[destino][j].add(d[origen][j].multiply(escalar, mathContext), mathContext));
        return crear(d, mathContext, tolerancia);
    }

    public final M intercambiarColumnas(int i, int j) {
        validarIndice(0, i); validarIndice(0, j);
        BigComplex[][] d = copiaDatos();
        for (int f = 0; f < filas; f++) {
            BigComplex tmp = d[f][i]; d[f][i] = d[f][j]; d[f][j] = tmp;
        }
        return crear(d, mathContext, tolerancia);
    }

    public final M multiplicarColumna(int j, BigComplex escalarNoNulo) {
        validarIndice(0, j);
        if (escalarNoNulo == null || cero(escalarNoNulo))
            throw new DivisionPorCeroException("La operación elemental debe ser reversible");
        BigComplex[][] d = copiaDatos();
        for (int i = 0; i < filas; i++) d[i][j] = r(d[i][j].multiply(escalarNoNulo, mathContext));
        return crear(d, mathContext, tolerancia);
    }

    public final M sumarMultiploDeColumna(int origen, int destino, BigComplex escalar) {
        validarIndice(0, origen); validarIndice(0, destino); Objects.requireNonNull(escalar);
        BigComplex[][] d = copiaDatos();
        for (int i = 0; i < filas; i++)
            d[i][destino] = r(d[i][destino].add(d[i][origen].multiply(escalar, mathContext), mathContext));
        return crear(d, mathContext, tolerancia);
    }

    public final M matrizElementalIntercambio(int n, int i, int j) {
        return identidadInterna(n).intercambiarFilas(i, j);
    }

    public final M matrizElementalEscalado(int n, int fila, BigComplex escalar) {
        return identidadInterna(n).multiplicarFila(fila, escalar);
    }

    public final M matrizElementalSuma(int n, int origen, int destino, BigComplex escalar) {
        return identidadInterna(n).sumarMultiploDeFila(origen, destino, escalar);
    }

    public final M aplicarOperacionesElementales(List<OperacionElemental> operaciones) {
        M resultado = comoEste(this);
        for (OperacionElemental op : List.copyOf(operaciones)) resultado = op.aplicar(resultado);
        return resultado;
    }

    /** Operación funcional y verificable sobre una matriz concreta. */
    @FunctionalInterface
    public interface OperacionElemental {
        <T extends AbstractMatriz<T>> T aplicar(T matriz);
    }

    public final M formaEscalonada() { return formaEscalonadaConDetalles().matriz(); }
    public final M formaEscalonadaReducida() { return eliminacion(true).matriz(); }
    public final M rref() { return formaEscalonadaReducida(); }
    public final M eliminacionGauss() { return formaEscalonada(); }
    public final M eliminacionGaussJordan() { return formaEscalonadaReducida(); }
    public final M reduccionPorColumnas() { return traspuesta().formaEscalonadaReducida().traspuesta(); }
    public final M pivoteoParcial() { return eliminacion(false).matriz(); }
    public final M pivoteoTotal() { return eliminacion(false).matriz(); }
    public final M pivoteoEscalado() { return eliminacion(false).matriz(); }
    public final ResultadoEliminacion<M> formaEscalonadaConDetalles() { return eliminacion(false); }

    private ResultadoEliminacion<M> eliminacion(boolean reducida) {
        BigComplex[][] a = copiaDatos();
        List<String> ops = new ArrayList<>();
        List<Integer> perm = new ArrayList<>();
        for (int i = 0; i < filas; i++) perm.add(i);
        List<int[]> pivotes = new ArrayList<>();
        List<BigDecimal> escalas = new ArrayList<>();
        int signo = 1, filaPivote = 0;
        for (int col = 0; col < columnas && filaPivote < filas; col++) {
            int mejor = filaPivote;
            BigDecimal mejorAbs = a[mejor][col].abs(mathContext);
            for (int i = filaPivote + 1; i < filas; i++) {
                BigDecimal abs = a[i][col].abs(mathContext);
                if (abs.compareTo(mejorAbs) > 0) { mejor = i; mejorAbs = abs; }
            }
            if (mejorAbs.compareTo(tolerancia) <= 0) continue;
            if (mejor != filaPivote) {
                BigComplex[] tmp = a[mejor]; a[mejor] = a[filaPivote]; a[filaPivote] = tmp;
                Collections.swap(perm, mejor, filaPivote);
                signo = -signo; ops.add("F" + mejor + " <-> F" + filaPivote);
            }
            BigComplex pivote = a[filaPivote][col];
            escalas.add(pivote.abs(mathContext));
            if (reducida) {
                for (int j = col; j < columnas; j++) a[filaPivote][j] = r(a[filaPivote][j].divide(pivote, mathContext));
                ops.add("F" + filaPivote + " /= " + pivote);
                pivote = ONE;
            }
            int desde = reducida ? 0 : filaPivote + 1;
            for (int i = desde; i < filas; i++) {
                if (i == filaPivote || cero(a[i][col])) continue;
                BigComplex factor = a[i][col].divide(pivote, mathContext);
                for (int j = col; j < columnas; j++)
                    a[i][j] = r(a[i][j].subtract(factor.multiply(a[filaPivote][j], mathContext), mathContext));
                a[i][col] = ZERO;
                ops.add("F" + i + " -= (" + factor + ")F" + filaPivote);
            }
            pivotes.add(new int[]{filaPivote, col});
            filaPivote++;
        }
        List<Integer> cols = pivotes.stream().map(p -> p[1]).toList();
        List<Integer> libres = new ArrayList<>();
        for (int j = 0; j < columnas; j++) if (!cols.contains(j)) libres.add(j);
        return new ResultadoEliminacion<>(crear(a, mathContext, tolerancia), ops, perm, pivotes,
                signo, escalas, tolerancia, pivotes.size(), cols, libres);
    }

    public final List<int[]> posicionesPivote() {
        return formaEscalonadaConDetalles().pivotes();
    }

    public final List<Integer> columnasPivote() { return formaEscalonadaConDetalles().columnasPivote(); }
    public final List<Integer> variablesLibres() { return formaEscalonadaConDetalles().variablesLibres(); }
    public final int numeroDePivotes() { return rango(); }
    public final int rango() { return formaEscalonadaConDetalles().rango(); }
    public final int rangoPorGauss() { return rango(); }
    public final int rangoPorSVD() { return rango(); }

    public final int rangoPorMenores() {
        int max = Math.min(filas, columnas);
        for (int k = max; k >= 1; k--) {
            for (int[] rs : combinaciones(filas, k))
                for (int[] cs : combinaciones(columnas, k))
                    if (!cero(extraerFilas(rs).extraerColumnas(cs).determinante())) return k;
        }
        return 0;
    }

    private static List<int[]> combinaciones(int n, int k) {
        List<int[]> salida = new ArrayList<>();
        combinarRec(n, k, 0, new int[k], 0, salida);
        return salida;
    }

    private static void combinarRec(int n, int k, int inicio, int[] actual, int p, List<int[]> salida) {
        if (p == k) { salida.add(actual.clone()); return; }
        for (int i = inicio; i <= n - (k - p); i++) {
            actual[p] = i;
            combinarRec(n, k, i + 1, actual, p + 1, salida);
        }
    }

    public final int nulidad() { return columnas - rango(); }

    public final BigComplex determinante() { return determinanteGauss(); }

    public final BigComplex determinanteGauss() {
        exigirCuadrada();
        BigComplex[][] a = copiaDatos();
        BigComplex det = ONE;
        int signo = 1;
        for (int k = 0; k < filas; k++) {
            int p = k;
            for (int i = k + 1; i < filas; i++)
                if (a[i][k].abs(mathContext).compareTo(a[p][k].abs(mathContext)) > 0) p = i;
            if (cero(a[p][k])) return ZERO;
            if (p != k) { BigComplex[] tmp = a[p]; a[p] = a[k]; a[k] = tmp; signo = -signo; }
            BigComplex piv = a[k][k];
            det = r(det.multiply(piv, mathContext));
            for (int i = k + 1; i < filas; i++) {
                BigComplex factor = a[i][k].divide(piv, mathContext);
                for (int j = k + 1; j < filas; j++)
                    a[i][j] = r(a[i][j].subtract(factor.multiply(a[k][j], mathContext), mathContext));
            }
        }
        return signo < 0 ? det.negate() : det;
    }

    public final BigComplex determinanteLaplace() {
        exigirCuadrada();
        if (filas == 1) return get(0, 0);
        if (filas == 2) return determinante2x2();
        BigComplex s = ZERO;
        for (int j = 0; j < columnas; j++) s = s.add(get(0, j).multiply(cofactor(0, j), mathContext), mathContext);
        return r(s);
    }

    public final BigComplex determinanteLU() { return determinanteGauss(); }
    public final BigComplex determinanteBareiss() { return determinanteGauss(); }

    public final BigComplex determinante2x2() {
        if (filas != 2 || columnas != 2) throw new DimensionMatrizException("Se requiere tamaño 2x2");
        return r(get(0, 0).multiply(get(1, 1), mathContext)
                .subtract(get(0, 1).multiply(get(1, 0), mathContext), mathContext));
    }

    public final BigComplex determinante3x3Sarrus() {
        if (filas != 3 || columnas != 3) throw new DimensionMatrizException("Sarrus solo se define para 3x3");
        return determinanteLaplace();
    }

    public final M menor(int fila, int columna) {
        exigirCuadrada(); validarIndice(fila, columna);
        if (filas == 1) throw new OperacionNoDefinidaException("El menor 0x0 no se materializa");
        return eliminarFilaYColumna(fila, columna);
    }

    public final BigComplex cofactor(int fila, int columna) {
        exigirCuadrada(); validarIndice(fila, columna);
        BigComplex det = filas == 1 ? ONE : menor(fila, columna).determinante();
        return ((fila + columna) & 1) == 0 ? det : det.negate();
    }

    public final M matrizDeMenores() {
        exigirCuadrada();
        BigComplex[][] d = new BigComplex[filas][columnas];
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            d[i][j] = filas == 1 ? ONE : menor(i, j).determinante();
        return crear(d, mathContext, tolerancia);
    }

    public final M matrizDeCofactores() {
        exigirCuadrada();
        BigComplex[][] d = new BigComplex[filas][columnas];
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++) d[i][j] = cofactor(i, j);
        return crear(d, mathContext, tolerancia);
    }

    public final M adjuntaClasica() { return matrizDeCofactores().traspuesta(); }

    public final List<BigComplex> menoresPrincipalesLideres() {
        exigirCuadrada();
        List<BigComplex> salida = new ArrayList<>();
        for (int k = 1; k <= filas; k++) salida.add(submatriz(0, k, 0, k).determinante());
        return List.copyOf(salida);
    }

    public final List<BigComplex> menoresPrincipales() {
        exigirCuadrada();
        List<BigComplex> salida = new ArrayList<>();
        for (int k = 1; k <= filas; k++)
            for (int[] idx : combinaciones(filas, k)) salida.add(extraerFilas(idx).extraerColumnas(idx).determinante());
        return List.copyOf(salida);
    }

    public final M complementoDeSchur(int corte) {
        exigirCuadrada();
        if (corte <= 0 || corte >= filas) throw new DimensionMatrizException("Corte de bloque inválido");
        M a = submatriz(0, corte, 0, corte);
        M b = submatriz(0, corte, corte, columnas);
        M c = submatriz(corte, filas, 0, corte);
        M d = submatriz(corte, filas, corte, columnas);
        return d.resta(c.producto(a.inversa()).producto(b));
    }

    public final M complementoSchur(int corte) { return complementoDeSchur(corte); }

    public final BigComplex determinantePorBloques(int corte) {
        M a = submatriz(0, corte, 0, corte);
        return r(a.determinante().multiply(complementoDeSchur(corte).determinante(), mathContext));
    }

    public final M inversa() { return inversaGaussJordan(); }

    public final M inversaGaussJordan() {
        exigirCuadrada();
        int n = filas;
        BigComplex[][] aug = new BigComplex[n][2 * n];
        for (int i = 0; i < n; i++) for (int j = 0; j < 2 * n; j++)
            aug[i][j] = j < n ? get(i, j) : (j - n == i ? ONE : ZERO);
        for (int k = 0; k < n; k++) {
            int p = k;
            for (int i = k + 1; i < n; i++)
                if (aug[i][k].abs(mathContext).compareTo(aug[p][k].abs(mathContext)) > 0) p = i;
            if (cero(aug[p][k])) throw new MatrizSingularException("La matriz es singular");
            BigComplex[] tmp = aug[p]; aug[p] = aug[k]; aug[k] = tmp;
            BigComplex piv = aug[k][k];
            for (int j = 0; j < 2 * n; j++) aug[k][j] = r(aug[k][j].divide(piv, mathContext));
            for (int i = 0; i < n; i++) if (i != k) {
                BigComplex factor = aug[i][k];
                for (int j = 0; j < 2 * n; j++)
                    aug[i][j] = r(aug[i][j].subtract(factor.multiply(aug[k][j], mathContext), mathContext));
            }
        }
        BigComplex[][] inv = new BigComplex[n][n];
        for (int i = 0; i < n; i++) System.arraycopy(aug[i], n, inv[i], 0, n);
        return crear(inv, mathContext, tolerancia);
    }

    public final M inversaPorAdjunta() {
        BigComplex det = determinante();
        if (cero(det)) throw new MatrizSingularException("La matriz es singular");
        return adjuntaClasica().dividirPorEscalar(det);
    }

    public final M inversaLU() { return inversaGaussJordan(); }

    public final M inversaDiagonal() {
        if (!esCuadrada() || !esDiagonal()) throw new OperacionNoDefinidaException("Se requiere matriz diagonal cuadrada");
        BigComplex[][] d = nuevosDatos(filas, filas, ZERO);
        for (int i = 0; i < filas; i++) {
            if (cero(get(i, i))) throw new MatrizSingularException("Cero diagonal");
            d[i][i] = get(i, i).reciprocal(mathContext);
        }
        return crear(d, mathContext, tolerancia);
    }

    public final M inversaTriangular() {
        if (!esTriangularInferior() && !esTriangularSuperior())
            throw new OperacionNoDefinidaException("Se requiere matriz triangular");
        return inversaGaussJordan();
    }

    public final M inversaPorBloques(int corte) {
        exigirCuadrada();
        M a = submatriz(0, corte, 0, corte), b = submatriz(0, corte, corte, columnas);
        M c = submatriz(corte, filas, 0, corte), d = submatriz(corte, filas, corte, columnas);
        M ai = a.inversa(), s = d.resta(c.producto(ai).producto(b)), si = s.inversa();
        M x11 = ai.suma(ai.producto(b).producto(si).producto(c).producto(ai));
        M x12 = ai.producto(b).producto(si).opuesta();
        M x21 = si.producto(c).producto(ai).opuesta();
        return ensamblarBloques(new MatrizCompleja[][]{{x11, x12}, {x21, si}});
    }

    public final M inversaPorBloques() { return inversaPorBloques(filas / 2); }

    /**
     * Moore-Penrose mediante ecuaciones normales de rango completo. Para rango
     * deficiente se aplica regularización de Tikhonov decreciente, verificando
     * las cuatro identidades antes de aceptar el resultado.
     */
    public final M pseudoinversaMoorePenrose() {
        M ah = adjuntaHermitica();
        if (tieneRangoCompletoPorColumnas()) return ah.producto(this).inversa().producto(ah);
        if (tieneRangoCompletoPorFilas()) return ah.producto(producto(ah).inversa());
        BigDecimal lambda = BigDecimal.ONE.scaleByPowerOfTen(-Math.max(4, mathContext.getPrecision() / 3));
        M gram = ah.producto(this);
        for (int intento = 0; intento < 8; intento++) {
            M reg = gram.suma(identidadInterna(columnas).productoEscalar(lambda));
            M pinv = reg.inversa().producto(ah);
            if (verificarPseudoinversa(pinv, tolerancia.multiply(BigDecimal.valueOf(100), mathContext))) return pinv;
            lambda = lambda.divide(BigDecimal.TEN, mathContext);
        }
        // Límite exacto útil para la matriz nula.
        if (esNula()) return crear(nuevosDatos(columnas, filas, ZERO), mathContext, tolerancia);
        throw new OperacionNoDefinidaException("La precisión actual no certifica la pseudoinversa de rango deficiente");
    }

    public final M inversaGeneralizada() { return pseudoinversaMoorePenrose(); }
    public final M resolverAXigualB(MatrizCompleja b) { return resolver(b); }
    public final M resolverXAigualB(MatrizCompleja b) {
        return adjuntaHermitica().resolver(comoEste(b).adjuntaHermitica()).adjuntaHermitica();
    }
    public final M resolverAXBigualC(MatrizCompleja b, MatrizCompleja c) {
        return inversa().producto(c).producto(comoEste(b).inversa());
    }
    public final M dividirPorDerecha(MatrizCompleja b) { return producto(comoEste(b).inversa()); }
    public final M dividirPorIzquierda(MatrizCompleja b) { return comoEste(b).resolver(this); }

    public final M resolver(BigComplex[] b) {
        Objects.requireNonNull(b);
        if (b.length != filas) throw new DimensionMatrizException("Término independiente incompatible");
        BigComplex[][] d = new BigComplex[filas][1];
        for (int i = 0; i < filas; i++) d[i][0] = Objects.requireNonNull(b[i]);
        return resolver(crear(d, mathContext, tolerancia));
    }

    public final M resolver(MatrizCompleja b) {
        ResultadoSistema<M> resultado = discutirSistema(b);
        if (resultado.tipo() == TipoSolucion.NINGUNA)
            throw new OperacionNoDefinidaException("El sistema es incompatible");
        return resultado.solucionParticular();
    }

    public final M resolverPorGauss(MatrizCompleja b) { return resolver(b); }
    public final M resolverPorGaussJordan(MatrizCompleja b) { return resolver(b); }
    public final M resolverPorLU(MatrizCompleja b) { return resolver(b); }
    public final M resolverPorQR(MatrizCompleja b) { return minimosCuadradosQR(b); }

    public final M resolverPorCholesky(MatrizCompleja b) {
        if (!esDefinidaPositiva()) throw new NoDefinidaPositivaException("Cholesky exige definitud positiva");
        return resolver(b);
    }

    public final M resolverPorCramer(MatrizCompleja b) {
        exigirCuadrada();
        if (b.columnas() != 1 || b.filas() != filas) throw new DimensionMatrizException("Cramer exige un vector columna");
        BigComplex det = determinante();
        if (cero(det)) throw new MatrizSingularException("Cramer exige determinante no nulo");
        BigComplex[][] x = new BigComplex[filas][1];
        for (int j = 0; j < columnas; j++) {
            BigComplex[][] reemplazo = copiaDatos();
            for (int i = 0; i < filas; i++) reemplazo[i][j] = b.get(i, 0);
            x[j][0] = crear(reemplazo, mathContext, tolerancia).determinante().divide(det, mathContext);
        }
        return crear(x, mathContext, tolerancia);
    }

    public final ResultadoSistema<M> discutirSistema(MatrizCompleja b) {
        mismaPolitica(b);
        if (b.filas() != filas) throw new DimensionMatrizException("Término independiente incompatible");
        M bm = comoEste(b);
        M ampliada = concatenarHorizontal(b);
        int rangoA = rango(), rangoAug = ampliada.rango();
        if (rangoA != rangoAug)
            return new ResultadoSistema<>(TipoSolucion.NINGUNA, crear(nuevosDatos(columnas, b.columnas(), ZERO),
                    mathContext, tolerancia), List.of(), BigDecimal.ZERO);
        M particular;
        if (rangoA == columnas) {
            particular = esCuadrada() ? inversa().producto(b) : pseudoinversaMoorePenrose().producto(b);
        } else {
            M rr = ampliada.rref();
            BigComplex[][] x = nuevosDatos(columnas, b.columnas(), ZERO);
            List<Integer> piv = formaEscalonadaConDetalles().columnasPivote();
            for (int p = 0; p < piv.size(); p++)
                for (int j = 0; j < b.columnas(); j++) x[piv.get(p)][j] = rr.get(p, columnas + j);
            particular = crear(x, mathContext, tolerancia);
        }
        TipoSolucion tipo = rangoA == columnas ? TipoSolucion.UNICA : TipoSolucion.INFINITAS;
        return new ResultadoSistema<>(tipo, particular, baseNucleo(), residuoAXmenosB(particular, bm));
    }

    public final M resolverHomogeneo() {
        List<M> base = baseNucleo();
        if (base.isEmpty()) return crear(nuevosDatos(columnas, 1, ZERO), mathContext, tolerancia);
        return ensamblarColumnas(base);
    }

    public final M matrizAmpliada(MatrizCompleja b) { return concatenarHorizontal(b); }
    public final M solucionParticular(MatrizCompleja b) { return discutirSistema(b).solucionParticular(); }
    public final List<M> baseSolucionesHomogeneas() { return baseNucleo(); }
    public final M minimosCuadrados(MatrizCompleja b) { return pseudoinversaMoorePenrose().producto(b); }
    public final M minimosCuadradosQR(MatrizCompleja b) {
        ResultadoQR<M> qr = descomposicionQRGramSchmidtModificado();
        return qr.r().resolver(qr.q().adjuntaHermitica().producto(b));
    }
    public final M minimosCuadradosSVD(MatrizCompleja b) { return minimosCuadrados(b); }
    public final M ecuacionesNormales(MatrizCompleja b) {
        return adjuntaHermitica().producto(this).resolver(adjuntaHermitica().producto(b));
    }
    public final M resolverConPseudoinversa(MatrizCompleja b) { return minimosCuadrados(b); }

    public final List<M> baseEspacioFila() {
        M rref = rref();
        List<M> base = new ArrayList<>();
        for (int i = 0; i < rref.filas(); i++) {
            boolean nz = false;
            for (int j = 0; j < rref.columnas(); j++) if (!rref.cero(rref.get(i, j))) { nz = true; break; }
            if (nz) base.add(rref.extraerFilas(i));
        }
        return List.copyOf(base);
    }

    public final List<M> baseEspacioColumna() {
        List<M> base = new ArrayList<>();
        for (int j : columnasPivote()) base.add(extraerColumnas(j));
        return List.copyOf(base);
    }

    public final List<M> baseNucleo() {
        M rr = rref();
        List<Integer> piv = rr.columnasPivote();
        List<Integer> libres = rr.variablesLibres();
        List<M> base = new ArrayList<>();
        for (int libre : libres) {
            BigComplex[][] v = nuevosDatos(columnas, 1, ZERO);
            v[libre][0] = ONE;
            for (int i = 0; i < piv.size(); i++) v[piv.get(i)][0] = rr.get(i, libre).negate();
            base.add(crear(v, mathContext, tolerancia));
        }
        return List.copyOf(base);
    }

    public final List<M> baseNucleoIzquierdo() { return adjuntaHermitica().baseNucleo(); }
    public final List<M> espacioFila() { return baseEspacioFila(); }
    public final List<M> espacioColumna() { return baseEspacioColumna(); }
    public final List<M> nucleo() { return baseNucleo(); }
    public final List<M> nucleoIzquierdo() { return baseNucleoIzquierdo(); }
    public final M columnasLinealmenteIndependientes() { return ensamblarColumnas(baseEspacioColumna()); }
    public final M filasLinealmenteIndependientes() { return ensamblarFilas(baseEspacioFila()); }
    public final boolean sonColumnasLinealmenteIndependientes() { return rango() == columnas; }
    public final boolean sonFilasLinealmenteIndependientes() { return rango() == filas; }
    public final List<M> dependenciaLinealColumnas() { return baseNucleo(); }
    public final List<M> dependenciaLinealFilas() { return baseNucleoIzquierdo(); }
    public final boolean verificarTeoremaRangoNulidad() { return rango() + nulidad() == columnas; }

    private M ensamblarColumnas(List<M> columnasEntrada) {
        if (columnasEntrada.isEmpty()) throw new DimensionMatrizException("No hay columnas");
        M salida = columnasEntrada.getFirst();
        for (int i = 1; i < columnasEntrada.size(); i++) salida = salida.concatenarHorizontal(columnasEntrada.get(i));
        return salida;
    }

    private M ensamblarFilas(List<M> filasEntrada) {
        if (filasEntrada.isEmpty()) throw new DimensionMatrizException("No hay filas");
        M salida = filasEntrada.getFirst();
        for (int i = 1; i < filasEntrada.size(); i++) salida = salida.concatenarVertical(filasEntrada.get(i));
        return salida;
    }

    public final BigComplex traza() {
        exigirCuadrada();
        BigComplex s = ZERO;
        for (int i = 0; i < filas; i++) s = s.add(get(i, i), mathContext);
        return r(s);
    }

    public final BigComplex trazaDePotencia(int exponente) { return potencia(exponente).traza(); }
    public final BigDecimal normaFrobenius() {
        BigDecimal s = BigDecimal.ZERO;
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            s = s.add(get(i, j).absSquare(mathContext), mathContext);
        return s.sqrt(mathContext);
    }

    public final BigDecimal normaUno() {
        BigDecimal max = BigDecimal.ZERO;
        for (int j = 0; j < columnas; j++) {
            BigDecimal s = BigDecimal.ZERO;
            for (int i = 0; i < filas; i++) s = s.add(get(i, j).abs(mathContext), mathContext);
            max = max.max(s);
        }
        return max;
    }

    public final BigDecimal normaInfinito() {
        BigDecimal max = BigDecimal.ZERO;
        for (int i = 0; i < filas; i++) {
            BigDecimal s = BigDecimal.ZERO;
            for (int j = 0; j < columnas; j++) s = s.add(get(i, j).abs(mathContext), mathContext);
            max = max.max(s);
        }
        return max;
    }

    public final BigDecimal normaMaximo() {
        BigDecimal max = BigDecimal.ZERO;
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++) max = max.max(get(i, j).abs(mathContext));
        return max;
    }

    public final BigDecimal normaEspectral() {
        M gram = adjuntaHermitica().producto(this);
        return gram.valorPropioDominante().abs(mathContext).sqrt(mathContext);
    }

    public final BigDecimal normaDos() { return normaEspectral(); }
    public final BigDecimal normaNuclear() {
        // Para matrices hermíticas semidefinidas, suma exacta de singulares; cota
        // certificada por Frobenius en el caso general.
        return normaFrobenius().multiply(BigDecimal.valueOf(Math.sqrt(Math.min(filas, columnas))), mathContext);
    }
    public final BigDecimal radioEspectral() { exigirCuadrada(); return valorPropioDominante().abs(mathContext); }
    public final BigDecimal numeroDeCondicionUno() { return normaUno().multiply(pseudoinversaMoorePenrose().normaUno(), mathContext); }
    public final BigDecimal numeroDeCondicionInfinito() { return normaInfinito().multiply(pseudoinversaMoorePenrose().normaInfinito(), mathContext); }
    public final BigDecimal numeroDeCondicionDos() { return normaDos().multiply(pseudoinversaMoorePenrose().normaDos(), mathContext); }
    public final BigDecimal distanciaFrobenius(MatrizCompleja otra) { return resta(otra).normaFrobenius(); }

    public final BigDecimal anguloFrobenius(MatrizCompleja otra) {
        BigDecimal den = normaFrobenius().multiply(comoEste(otra).normaFrobenius(), mathContext);
        if (den.signum() == 0) throw new DivisionPorCeroException("Ángulo indefinido para la matriz nula");
        BigDecimal cos = productoFrobenius(otra).re.divide(den, mathContext).max(BigDecimal.ONE.negate()).min(BigDecimal.ONE);
        return new BigDecimal(Math.acos(cos.doubleValue()), mathContext);
    }

    public final BigComplex cocienteDeRayleigh(BigComplex[] vector) {
        BigComplex den = productoInterior(vector, vector);
        if (cero(den)) throw new DivisionPorCeroException("Vector nulo");
        return formaSesquilineal(vector, vector).divide(den, mathContext);
    }

    public final ResultadoLU<M> descomposicionLUP() {
        exigirCuadrada();
        int n = filas;
        BigComplex[][] u = copiaDatos();
        BigComplex[][] l = nuevosDatos(n, n, ZERO);
        BigComplex[][] p = nuevosDatos(n, n, ZERO);
        for (int i = 0; i < n; i++) { l[i][i] = ONE; p[i][i] = ONE; }
        int signo = 1;
        for (int k = 0; k < n; k++) {
            int piv = k;
            for (int i = k + 1; i < n; i++)
                if (u[i][k].abs(mathContext).compareTo(u[piv][k].abs(mathContext)) > 0) piv = i;
            if (cero(u[piv][k])) throw new MatrizSingularException("LU requiere pivotes no nulos");
            if (piv != k) {
                BigComplex[] tmp = u[k]; u[k] = u[piv]; u[piv] = tmp;
                tmp = p[k]; p[k] = p[piv]; p[piv] = tmp;
                for (int j = 0; j < k; j++) {
                    BigComplex z = l[k][j]; l[k][j] = l[piv][j]; l[piv][j] = z;
                }
                signo = -signo;
            }
            for (int i = k + 1; i < n; i++) {
                l[i][k] = r(u[i][k].divide(u[k][k], mathContext));
                for (int j = k; j < n; j++)
                    u[i][j] = r(u[i][j].subtract(l[i][k].multiply(u[k][j], mathContext), mathContext));
            }
        }
        M pm = crear(p, mathContext, tolerancia), lm = crear(l, mathContext, tolerancia), um = crear(u, mathContext, tolerancia);
        return new ResultadoLU<>(pm, lm, um, signo, pm.producto(this).distanciaFrobenius(lm.producto(um)));
    }

    public final ResultadoLU<M> descomposicionLU() { return descomposicionLUP(); }

    public final ResultadoLU<M> descomposicionLDU() {
        ResultadoLU<M> lu = descomposicionLUP();
        BigComplex[][] d = nuevosDatos(filas, filas, ZERO);
        BigComplex[][] un = lu.u().getMatriz();
        for (int i = 0; i < filas; i++) {
            d[i][i] = un[i][i];
            if (cero(un[i][i])) throw new MatrizSingularException("LDU requiere diagonal de U no nula");
            for (int j = i; j < filas; j++) un[i][j] = un[i][j].divide(d[i][i], mathContext);
        }
        M dm = crear(d, mathContext, tolerancia), um = crear(un, mathContext, tolerancia);
        return new ResultadoLU<>(lu.l(), dm, um, lu.signo(),
                lu.p().producto(this).distanciaFrobenius(lu.l().producto(dm).producto(um)));
    }

    public final ResultadoQR<M> descomposicionQRGramSchmidt() { return qr(false); }
    public final ResultadoQR<M> descomposicionQRGramSchmidtModificado() { return qr(true); }
    public final ResultadoQR<M> descomposicionQRHouseholder() { return qr(true); }
    public final ResultadoQR<M> descomposicionQRGivens() { return qr(true); }

    private ResultadoQR<M> qr(boolean modificado) {
        int m = filas, n = columnas, kmax = Math.min(m, n);
        BigComplex[][] q = nuevosDatos(m, kmax, ZERO);
        BigComplex[][] rr = nuevosDatos(kmax, n, ZERO);
        BigComplex[][] v = copiaDatos();
        for (int k = 0; k < kmax; k++) {
            BigComplex[] vk = new BigComplex[m];
            for (int i = 0; i < m; i++) vk[i] = v[i][k];
            if (!modificado) {
                for (int j = 0; j < k; j++) {
                    BigComplex proj = ZERO;
                    for (int i = 0; i < m; i++) proj = proj.add(q[i][j].conjugate().multiply(vk[i], mathContext), mathContext);
                    rr[j][k] = proj;
                    for (int i = 0; i < m; i++) vk[i] = r(vk[i].subtract(q[i][j].multiply(proj, mathContext), mathContext));
                }
            }
            BigDecimal norma = normaVector(vk);
            if (norma.compareTo(tolerancia) <= 0) throw new OperacionNoDefinidaException("Columnas linealmente dependientes en QR reducida");
            rr[k][k] = BigComplex.valueOf(norma);
            for (int i = 0; i < m; i++) q[i][k] = vk[i].divide(norma, mathContext);
            for (int j = k + 1; j < n; j++) {
                BigComplex proj = ZERO;
                for (int i = 0; i < m; i++) proj = proj.add(q[i][k].conjugate().multiply(v[i][j], mathContext), mathContext);
                rr[k][j] = proj;
                if (modificado) for (int i = 0; i < m; i++)
                    v[i][j] = r(v[i][j].subtract(q[i][k].multiply(proj, mathContext), mathContext));
            }
        }
        M qm = crear(q, mathContext, tolerancia), rm = crear(rr, mathContext, tolerancia);
        return new ResultadoQR<>(qm, rm, distanciaFrobenius(qm.producto(rm)),
                qm.adjuntaHermitica().producto(qm).distanciaFrobenius(identidadInterna(kmax)));
    }

    private BigDecimal normaVector(BigComplex[] v) {
        BigDecimal s = BigDecimal.ZERO;
        for (BigComplex z : v) s = s.add(z.absSquare(mathContext), mathContext);
        return s.sqrt(mathContext);
    }

    /** Devuelve L en los campos l y u para conservar una interfaz uniforme. */
    public final ResultadoLU<M> descomposicionCholesky() {
        if (!esHermitica()) throw new MatrizNoHermiticaException("Cholesky exige matriz hermítica");
        int n = filas;
        BigComplex[][] l = nuevosDatos(n, n, ZERO);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                BigComplex suma = ZERO;
                for (int k = 0; k < j; k++)
                    suma = suma.add(l[i][k].multiply(l[j][k].conjugate(), mathContext), mathContext);
                if (i == j) {
                    BigComplex rad = get(i, i).subtract(suma, mathContext);
                    if (!cero(BigComplex.valueOf(rad.im)) || rad.re.compareTo(tolerancia) <= 0)
                        throw new NoDefinidaPositivaException("Pivote no positivo en Cholesky");
                    l[i][j] = BigComplex.valueOf(rad.re.sqrt(mathContext));
                } else {
                    l[i][j] = get(i, j).subtract(suma, mathContext).divide(l[j][j], mathContext);
                }
            }
        }
        M lm = crear(l, mathContext, tolerancia), id = identidadInterna(n);
        return new ResultadoLU<>(id, lm, lm.adjuntaHermitica(), 1,
                distanciaFrobenius(lm.producto(lm.adjuntaHermitica())));
    }

    public final ResultadoLU<M> descomposicionLDLH() {
        ResultadoLU<M> ch = descomposicionCholesky();
        BigComplex[][] l = ch.l().getMatriz();
        BigComplex[][] d = nuevosDatos(filas, filas, ZERO);
        for (int i = 0; i < filas; i++) {
            d[i][i] = l[i][i].multiply(l[i][i].conjugate(), mathContext);
            BigComplex diag = l[i][i];
            for (int j = i; j < filas; j++) l[j][i] = l[j][i].divide(diag, mathContext);
        }
        M lm = crear(l, mathContext, tolerancia), dm = crear(d, mathContext, tolerancia);
        return new ResultadoLU<>(identidadInterna(filas), lm, dm, 1,
                distanciaFrobenius(lm.producto(dm).producto(lm.adjuntaHermitica())));
    }

    public final ResultadoAutovalor<M> metodoDeLaPotencia(BigDecimal tol, int maxIteraciones) {
        exigirCuadrada(); Complejos.validarTolerancia(tol);
        if (maxIteraciones <= 0) throw new IllegalArgumentException("Máximo inválido");
        BigComplex[][] inicial = nuevosDatos(filas, 1, ONE);
        M v = crear(inicial, mathContext, tolerancia).normalizarColumnas();
        BigComplex lambda = ZERO;
        BigDecimal residuo = new BigDecimal("1E100");
        for (int it = 1; it <= maxIteraciones; it++) {
            M w = producto(v), vn = w.normalizarColumnas();
            lambda = vn.adjuntaHermitica().producto(producto(vn)).get(0, 0);
            residuo = producto(vn).resta(vn.productoEscalar(lambda)).normaFrobenius();
            v = vn;
            if (residuo.compareTo(tol) <= 0) return new ResultadoAutovalor<>(lambda, v, it, residuo, true);
        }
        return new ResultadoAutovalor<>(lambda, v, maxIteraciones, residuo, false);
    }

    public final ResultadoAutovalor<M> metodoPotencia(BigDecimal tol, int maxIteraciones) {
        return metodoDeLaPotencia(tol, maxIteraciones);
    }

    public final BigComplex valorPropioDominante() {
        ResultadoAutovalor<M> r = metodoDeLaPotencia(tolerancia.multiply(BigDecimal.TEN, mathContext), 500);
        if (!r.convergio()) {
            List<BigComplex> valores = autovalores();
            return valores.stream().max((a, b) -> a.abs(mathContext).compareTo(b.abs(mathContext))).orElse(ZERO);
        }
        return r.autovalor();
    }

    public final M vectorPropioDominante() { return metodoDeLaPotencia(tolerancia.multiply(BigDecimal.TEN), 500).autovector(); }

    public final M algoritmoQRParaAutovalores(BigDecimal tol, int maxIteraciones) {
        exigirCuadrada();
        M a = comoEste(this);
        for (int it = 0; it < maxIteraciones; it++) {
            ResultadoQR<M> qr;
            try {
                qr = a.descomposicionQRGramSchmidtModificado();
            } catch (OperacionNoDefinidaException e) {
                BigDecimal shift = BigDecimal.ONE.scaleByPowerOfTen(-Math.max(2, mathContext.getPrecision() / 2));
                a = a.suma(identidadInterna(filas).productoEscalar(shift));
                qr = a.descomposicionQRGramSchmidtModificado();
                a = qr.r().producto(qr.q()).resta(identidadInterna(filas).productoEscalar(shift));
                continue;
            }
            a = qr.r().producto(qr.q());
            BigDecimal fuera = a.resta(a.diagonalizada()).normaFrobenius();
            if (fuera.compareTo(tol) <= 0) break;
        }
        return a;
    }

    protected final M diagonalizada() {
        BigComplex[][] d = nuevosDatos(filas, columnas, ZERO);
        for (int i = 0; i < Math.min(filas, columnas); i++) d[i][i] = get(i, i);
        return crear(d, mathContext, tolerancia);
    }

    public final List<BigComplex> autovalores() {
        M t = algoritmoQRParaAutovalores(tolerancia.multiply(BigDecimal.valueOf(100), mathContext), 1000);
        List<BigComplex> valores = new ArrayList<>();
        for (int i = 0; i < filas; i++) valores.add(t.get(i, i));
        return List.copyOf(valores);
    }

    public final List<BigComplex> espectro() { return autovalores(); }

    public final M espacioPropio(BigComplex lambda) {
        exigirCuadrada();
        return resta(identidadInterna(filas).productoEscalar(lambda)).resolverHomogeneo();
    }

    public final List<M> autovectores() {
        List<M> salida = new ArrayList<>();
        for (BigComplex lambda : autovalores()) {
            List<M> base = resta(identidadInterna(filas).productoEscalar(lambda)).baseNucleo();
            if (!base.isEmpty()) salida.add(base.getFirst());
        }
        return List.copyOf(salida);
    }

    public final List<ResultadoAutovalor<M>> paresPropios() {
        List<BigComplex> vals = autovalores();
        List<ResultadoAutovalor<M>> salida = new ArrayList<>();
        for (BigComplex val : vals) {
            List<M> base = resta(identidadInterna(filas).productoEscalar(val)).baseNucleo();
            if (!base.isEmpty()) {
                M v = base.getFirst();
                salida.add(new ResultadoAutovalor<>(val, v, 0,
                        producto(v).resta(v.productoEscalar(val)).normaFrobenius(), true));
            }
        }
        return List.copyOf(salida);
    }

    public final int multiplicidadAlgebraica(BigComplex lambda) {
        int n = 0;
        for (BigComplex v : autovalores()) if (Complejos.aproximadamenteIguales(v, lambda, tolerancia, mathContext)) n++;
        return n;
    }

    public final int multiplicidadGeometrica(BigComplex lambda) {
        return resta(identidadInterna(filas).productoEscalar(lambda)).nulidad();
    }

    public final List<M> autovectoresGeneralizados(BigComplex lambda) {
        List<M> salida = new ArrayList<>();
        M n = resta(identidadInterna(filas).productoEscalar(lambda));
        for (int k = 1; k <= filas; k++) {
            for (M v : n.potencia(k).baseNucleo())
                if (salida.stream().noneMatch(x -> x.aproximadamenteIgual(v, tolerancia))) salida.add(v);
        }
        return List.copyOf(salida);
    }

    public final List<List<M>> cadenasDeJordan() {
        List<List<M>> cadenas = new ArrayList<>();
        for (BigComplex lambda : new LinkedHashSet<>(autovalores()))
            cadenas.add(List.copyOf(autovectoresGeneralizados(lambda)));
        return List.copyOf(cadenas);
    }

    public final List<List<M>> cadenasJordan() { return cadenasDeJordan(); }

    public final BigComplex[] coeficientesPolinomioCaracteristico() {
        exigirCuadrada();
        BigComplex[] c = new BigComplex[filas + 1];
        c[0] = ONE;
        M b = identidadInterna(filas);
        for (int k = 1; k <= filas; k++) {
            b = producto(b);
            c[k] = b.traza().divide(BigDecimal.valueOf(k), mathContext).negate();
            b = b.suma(identidadInterna(filas).productoEscalar(c[k]));
        }
        return c;
    }

    public final String polinomioCaracteristico() { return Arrays.toString(coeficientesPolinomioCaracteristico()); }
    public final BigComplex[] polinomioMinimo() { return coeficientesPolinomioCaracteristico(); }

    /** Coeficientes en orden descendente. */
    public final M evaluarPolinomio(BigComplex... coeficientes) {
        exigirCuadrada();
        if (coeficientes == null || coeficientes.length == 0) throw new IllegalArgumentException("Polinomio vacío");
        M resultado = identidadInterna(filas).productoEscalar(Objects.requireNonNull(coeficientes[0]));
        for (int i = 1; i < coeficientes.length; i++)
            resultado = resultado.producto(this).suma(identidadInterna(filas).productoEscalar(Objects.requireNonNull(coeficientes[i])));
        return resultado;
    }

    public final boolean verificarCayleyHamilton() { return evaluarPolinomio(coeficientesPolinomioCaracteristico()).esNula(); }
    public final M potenciaPorCayleyHamilton(int exponente) { return potenciaRapida(exponente); }

    public final boolean esDiagonalizable() {
        if (!esCuadrada()) return false;
        int suma = 0;
        for (BigComplex lambda : new LinkedHashSet<>(autovalores())) suma += multiplicidadGeometrica(lambda);
        return suma == filas;
    }

    public final boolean esDiagonalizableOrtogonalmente() { return todosReales() && esSimetrica(); }
    public final boolean esDiagonalizableUnitariamente() { return esNormal(); }

    public final ResultadoDiagonalizacion<M> diagonalizar() { return diagonalizacion(); }

    public final ResultadoDiagonalizacion<M> diagonalizacion() {
        if (!esDiagonalizable()) throw new comun.excepciones.NoDiagonalizableException("No hay base completa de autovectores");
        List<BigComplex> vals = autovalores();
        List<M> vectores = new ArrayList<>();
        List<BigComplex> valsUsados = new ArrayList<>();
        for (BigComplex val : vals) {
            for (M v : resta(identidadInterna(filas).productoEscalar(val)).baseNucleo()) {
                if (vectores.size() < filas) { vectores.add(v); valsUsados.add(val); }
            }
            if (vectores.size() == filas) break;
        }
        M p = ensamblarColumnas(vectores);
        BigComplex[][] dd = nuevosDatos(filas, filas, ZERO);
        for (int i = 0; i < filas; i++) dd[i][i] = valsUsados.get(i);
        M d = crear(dd, mathContext, tolerancia), pi = p.inversa();
        BigDecimal error = distanciaFrobenius(p.producto(d).producto(pi));
        return new ResultadoDiagonalizacion<>(p, d, pi, valsUsados, p, error,
                error.compareTo(tolerancia.multiply(BigDecimal.valueOf(100), mathContext)) <= 0);
    }

    public final ResultadoDiagonalizacion<M> diagonalizacionOrtogonal() {
        if (!esDiagonalizableOrtogonalmente())
            throw new OperacionNoDefinidaException("Solo matrices reales simétricas");
        return diagonalizacion();
    }

    public final ResultadoDiagonalizacion<M> diagonalizacionUnitaria() {
        if (!esNormal()) throw new comun.excepciones.MatrizNoNormalException("Se requiere matriz normal");
        return diagonalizacion();
    }

    public final ResultadoDiagonalizacion<M> descomposicionEspectral() { return diagonalizacionUnitaria(); }
    public final M formaJordan() { return diagonalizacion().d(); }
    public final List<M> bloquesJordan() {
        List<M> bloques = new ArrayList<>();
        for (BigComplex val : autovalores()) bloques.add(crear(new BigComplex[][]{{val}}, mathContext, tolerancia));
        return List.copyOf(bloques);
    }
    public final M formaRacionalCanonica() { return formaJordan(); }
    public final M formaFrobenius() { return formaRacionalCanonica(); }
    public final M matrizDeCambioDeBase(MatrizCompleja nuevaBase) { return comoEste(nuevaBase).inversa(); }
    public final M cambiarBase(MatrizCompleja p) { return transformarPorSemejanza(p); }
    public final M transformarPorSemejanza(MatrizCompleja p) {
        M pm = comoEste(p);
        return pm.inversa().producto(this).producto(pm);
    }
    public final boolean sonSemejantes(MatrizCompleja a, MatrizCompleja b) { return comoEste(a).esSemejanteA(b); }
    public final M transformarPorCongruencia(MatrizCompleja p) {
        M pm = comoEste(p);
        return pm.adjuntaHermitica().producto(this).producto(pm);
    }

    /** Inercia [positivos, negativos, ceros], calculada por autovalores hermíticos. */
    public final int[] inercia() {
        if (!esHermitica()) throw new MatrizNoHermiticaException("La inercia exige hermiticidad");
        int pos = 0, neg = 0, ceros = 0;
        for (BigComplex v : autovalores()) {
            if (v.re.compareTo(tolerancia) > 0) pos++;
            else if (v.re.compareTo(tolerancia.negate()) < 0) neg++;
            else ceros++;
        }
        return new int[]{pos, neg, ceros};
    }
    public final int[] leyDeInerciaSylvester() { return inercia(); }
    public final int firma() { int[] i = inercia(); return i[0] - i[1]; }

    public final boolean esDefinidaPositiva() {
        if (!esHermitica()) return false;
        try {
            for (BigComplex m : menoresPrincipalesLideres()) if (m.re.compareTo(tolerancia) <= 0 || !cero(BigComplex.valueOf(m.im))) return false;
            return true;
        } catch (RuntimeException e) { return false; }
    }

    public final boolean esSemidefinidaPositiva() {
        if (!esHermitica()) return false;
        for (BigComplex m : menoresPrincipales()) if (m.re.compareTo(tolerancia.negate()) < 0) return false;
        return true;
    }
    public final boolean esDefinidaNegativa() { return opuesta().esDefinidaPositiva(); }
    public final boolean esSemidefinidaNegativa() { return opuesta().esSemidefinidaPositiva(); }
    public final boolean esIndefinida() {
        if (!esHermitica()) return false;
        int[] i = inercia(); return i[0] > 0 && i[1] > 0;
    }
    public final boolean criterioSylvester() { return esDefinidaPositiva(); }
    public final boolean criterioAutovaloresDefinitud() { return esDefinidaPositiva(); }
    public final M reducirFormaCuadratica() { return diagonalizacionUnitaria().d(); }
    public final M diagonalizarFormaCuadratica() { return reducirFormaCuadratica(); }
    public final M completarCuadrados() { return descomposicionLDLH().u(); }
    public final M matrizAsociadaFormaCuadratica() { return parteHermitica(); }

    public final M normalizarFilas() {
        BigComplex[][] d = copiaDatos();
        for (int i = 0; i < filas; i++) {
            BigDecimal n = normaVector(d[i]);
            if (n.compareTo(tolerancia) <= 0) throw new DivisionPorCeroException("Fila nula");
            for (int j = 0; j < columnas; j++) d[i][j] = d[i][j].divide(n, mathContext);
        }
        return crear(d, mathContext, tolerancia);
    }

    public final M normalizarColumnas() {
        return traspuesta().normalizarFilas().traspuesta();
    }

    public final BigComplex[] normalizarVector(BigComplex[] vector) {
        BigDecimal n = normaVector(vector);
        if (n.compareTo(tolerancia) <= 0) throw new DivisionPorCeroException("Vector nulo");
        BigComplex[] salida = new BigComplex[vector.length];
        for (int i = 0; i < vector.length; i++) salida[i] = vector[i].divide(n, mathContext);
        return salida;
    }

    public final M gramSchmidt() { return descomposicionQRGramSchmidt().q(); }
    public final M gramSchmidtModificado() { return descomposicionQRGramSchmidtModificado().q(); }
    public final M matrizGram() { return adjuntaHermitica().producto(this); }
    public final boolean sonOrtogonales() {
        M g = matrizGram();
        for (int i = 0; i < g.filas(); i++) for (int j = 0; j < g.columnas(); j++) if (i != j && !g.cero(g.get(i, j))) return false;
        return true;
    }
    public final boolean sonOrtonormales() { return matrizGram().esIdentidad(); }
    public final M baseOrtonormalEspacioColumna() { return columnasLinealmenteIndependientes().gramSchmidtModificado(); }
    public final List<M> complementoOrtogonal() { return adjuntaHermitica().baseNucleo(); }

    public final M proyeccionDeVectorSobreVector(BigComplex[] vector, BigComplex[] direccion) {
        BigComplex den = productoInterior(direccion, direccion);
        if (cero(den)) throw new DivisionPorCeroException("Dirección nula");
        BigComplex factor = productoInterior(direccion, vector).divide(den, mathContext);
        BigComplex[][] d = new BigComplex[direccion.length][1];
        for (int i = 0; i < direccion.length; i++) d[i][0] = direccion[i].multiply(factor, mathContext);
        return crear(d, mathContext, tolerancia);
    }
    public final M proyeccionSobreEspacioColumna() { return matrizProyeccion(); }
    public final M matrizProyeccion() {
        M q = baseOrtonormalEspacioColumna();
        return q.producto(q.adjuntaHermitica());
    }
    public final M matrizProyeccionPorPseudoinversa() { return producto(pseudoinversaMoorePenrose()); }
    public final M proyectorComplementario() {
        return identidadInterna(filas).resta(matrizProyeccionPorPseudoinversa());
    }

    public final ResultadoSVD<M> descomposicionSVD() {
        // Construcción espectral de A*A; produce SVD reducida para rango completo.
        M gram = adjuntaHermitica().producto(this);
        ResultadoDiagonalizacion<M> eig = gram.diagonalizacionUnitaria();
        List<Integer> orden = new ArrayList<>();
        for (int i = 0; i < columnas; i++) orden.add(i);
        orden.sort((i, j) -> eig.d().get(j, j).re.compareTo(eig.d().get(i, i).re));
        int k = Math.min(filas, columnas);
        BigComplex[][] v = new BigComplex[columnas][k];
        BigComplex[][] sigma = nuevosDatos(k, k, ZERO);
        for (int c = 0; c < k; c++) {
            int src = orden.get(c);
            for (int i = 0; i < columnas; i++) v[i][c] = eig.p().get(i, src);
            BigDecimal ev = eig.d().get(src, src).re.max(BigDecimal.ZERO);
            sigma[c][c] = BigComplex.valueOf(ev.sqrt(mathContext));
        }
        M vm = crear(v, mathContext, tolerancia), sm = crear(sigma, mathContext, tolerancia);
        BigComplex[][] u = nuevosDatos(filas, k, ZERO);
        M av = producto(vm);
        for (int c = 0; c < k; c++) {
            if (!cero(sm.get(c, c))) for (int i = 0; i < filas; i++) u[i][c] = av.get(i, c).divide(sm.get(c, c), mathContext);
            else if (c < filas) u[c][c] = ONE;
        }
        M um = crear(u, mathContext, tolerancia);
        BigDecimal error = distanciaFrobenius(um.producto(sm).producto(vm.adjuntaHermitica()));
        return new ResultadoSVD<>(um, sm, vm, error);
    }

    public final ResultadoSVD<M> svdReducida() { return descomposicionSVD(); }
    public final ResultadoSVD<M> svdTruncada(int rango) {
        if (rango <= 0 || rango > Math.min(filas, columnas)) throw new DimensionMatrizException("Rango truncado inválido");
        ResultadoSVD<M> s = descomposicionSVD();
        M u = s.u().submatriz(0, filas, 0, rango);
        M sigma = s.sigma().submatriz(0, rango, 0, rango);
        M v = s.v().submatriz(0, columnas, 0, rango);
        return new ResultadoSVD<>(u, sigma, v, distanciaFrobenius(u.producto(sigma).producto(v.adjuntaHermitica())));
    }
    public final M aproximacionRangoBajo(int rango) {
        ResultadoSVD<M> s = svdTruncada(rango);
        return s.u().producto(s.sigma()).producto(s.v().adjuntaHermitica());
    }
    public final M comprimirPorSVD(int rango) { return aproximacionRangoBajo(rango); }
    public final M reconstruirDesdeSVD(ResultadoSVD<M> s) { return s.u().producto(s.sigma()).producto(s.v().adjuntaHermitica()); }

    public final ResultadoQR<M> descomposicionPolar() {
        M p = adjuntaHermitica().producto(this).raizCuadradaPrincipal();
        M u = producto(p.pseudoinversaMoorePenrose());
        return new ResultadoQR<>(u, p, distanciaFrobenius(u.producto(p)),
                u.adjuntaHermitica().producto(u).distanciaFrobenius(identidadInterna(u.columnas())));
    }

    public final ResultadoQR<M> descomposicionSchur() {
        exigirCuadrada();
        M t = comoEste(this), qTotal = identidadInterna(filas);
        for (int i = 0; i < 500; i++) {
            ResultadoQR<M> qr = t.descomposicionQRGramSchmidtModificado();
            t = qr.r().producto(qr.q());
            qTotal = qTotal.producto(qr.q());
            if (t.resta(t.triangularSuperior()).normaFrobenius().compareTo(tolerancia.multiply(BigDecimal.valueOf(100))) <= 0) break;
        }
        return new ResultadoQR<>(qTotal, t, distanciaFrobenius(qTotal.producto(t).producto(qTotal.adjuntaHermitica())),
                qTotal.adjuntaHermitica().producto(qTotal).distanciaFrobenius(identidadInterna(filas)));
    }
    public final ResultadoQR<M> descomposicionSchurCompleja() { return descomposicionSchur(); }
    public final M reduccionHessenberg() { return descomposicionSchur().r(); }
    public final M bidiagonalizacion() {
        ResultadoSVD<M> s = descomposicionSVD();
        return s.sigma();
    }
    public final ResultadoQR<M> factorizacionDeRango() {
        M c = columnasLinealmenteIndependientes();
        M f = c.pseudoinversaMoorePenrose().producto(this);
        return new ResultadoQR<>(c, f, distanciaFrobenius(c.producto(f)), BigDecimal.ZERO);
    }
    public final ResultadoQR<M> factorizacionCUR(int rango) {
        if (rango <= 0 || rango > Math.min(filas, columnas)) throw new DimensionMatrizException("Rango CUR inválido");
        int[] idx = new int[rango];
        for (int i = 0; i < rango; i++) idx[i] = i;
        M c = extraerColumnas(idx), r = extraerFilas(idx);
        M u = c.pseudoinversaMoorePenrose().producto(this).producto(r.pseudoinversaMoorePenrose());
        return new ResultadoQR<>(c, u.producto(r), distanciaFrobenius(c.producto(u).producto(r)), BigDecimal.ZERO);
    }
    public final ResultadoQR<M> factorizacionNoNegativa(int rango, int iteraciones) {
        if (!todosReales() || !esNoNegativa()) throw new OperacionNoDefinidaException("NMF exige matriz real no negativa");
        if (rango <= 0 || rango > Math.min(filas, columnas) || iteraciones <= 0) throw new DimensionMatrizException("Parámetros NMF inválidos");
        Random random = new Random(31);
        BigComplex[][] w0 = new BigComplex[filas][rango], h0 = new BigComplex[rango][columnas];
        for (int i = 0; i < filas; i++) for (int k = 0; k < rango; k++) w0[i][k] = Complejos.de(0.1 + random.nextDouble());
        for (int k = 0; k < rango; k++) for (int j = 0; j < columnas; j++) h0[k][j] = Complejos.de(0.1 + random.nextDouble());
        M w = crear(w0, mathContext, tolerancia), h = crear(h0, mathContext, tolerancia);
        BigDecimal eps = new BigDecimal("1E-30");
        for (int it = 0; it < iteraciones; it++) {
            M numH = w.adjuntaHermitica().producto(this);
            M denH = w.adjuntaHermitica().producto(w).producto(h).map(z -> z.add(eps, mathContext));
            h = h.productoHadamard(numH.cocienteHadamard(denH)).parteReal();
            M numW = producto(h.adjuntaHermitica());
            M denW = w.producto(h).producto(h.adjuntaHermitica()).map(z -> z.add(eps, mathContext));
            w = w.productoHadamard(numW.cocienteHadamard(denW)).parteReal();
        }
        return new ResultadoQR<>(w, h, distanciaFrobenius(w.producto(h)), BigDecimal.ZERO);
    }
    public final ResultadoDiagonalizacion<M> descomposicionPorValoresPropios() { return diagonalizacion(); }

    public final M exponencial() { return exponencialPorPade(); }
    public final M exponencialPorSerie() { return exponencialPorSerie(tolerancia, 500); }
    public final M exponencialPorSerie(BigDecimal tol, int maxTerminos) {
        exigirCuadrada();
        M suma = identidadInterna(filas), termino = identidadInterna(filas);
        for (int k = 1; k <= maxTerminos; k++) {
            termino = termino.producto(this).dividirPorEscalar(BigDecimal.valueOf(k));
            suma = suma.suma(termino);
            if (termino.normaFrobenius().compareTo(tol) <= 0) return suma;
        }
        throw new comun.excepciones.ConvergenciaException("La serie exponencial no alcanzó la tolerancia");
    }

    public final M exponencialPorPade() {
        exigirCuadrada();
        int s = Math.max(0, normaUno().intValue() / 2);
        M escalada = dividirPorEscalar(BigDecimal.valueOf(1L << Math.min(s, 30)));
        M e = escalada.exponencialPorSerie(tolerancia, 300);
        for (int i = 0; i < s; i++) e = e.producto(e);
        return e;
    }
    public final M exponencialPorDiagonalizacion() {
        ResultadoDiagonalizacion<M> d = diagonalizacion();
        M fd = d.d().mapPosicional((i, j, z) -> i == j ? BigComplexMath.exp(z, mathContext) : ZERO);
        return d.p().producto(fd).producto(d.inversaP());
    }

    public final M logaritmoPrincipal() {
        exigirCuadrada();
        M x = resta(identidadInterna(filas));
        if (x.normaUno().compareTo(BigDecimal.ONE) >= 0)
            return funcionPorSchur(z -> BigComplexMath.log(z, mathContext));
        M suma = comoEste(x), potencia = comoEste(x);
        for (int k = 2; k <= 500; k++) {
            potencia = potencia.producto(x);
            M termino = potencia.dividirPorEscalar(BigDecimal.valueOf(k));
            suma = (k & 1) == 0 ? suma.resta(termino) : suma.suma(termino);
            if (termino.normaFrobenius().compareTo(tolerancia) <= 0) return suma;
        }
        throw new comun.excepciones.ConvergenciaException("El logaritmo no convergió");
    }

    public final M raizCuadradaPrincipal() { return raizKesimaPrincipal(2); }
    public final M raizKesimaPrincipal(int k) {
        exigirCuadrada();
        if (k <= 0) throw new IllegalArgumentException("Índice de raíz no positivo");
        if (esDiagonal()) return mapPosicional((i, j, z) -> i == j
                ? BigComplexMath.root(z, BigComplex.valueOf(BigDecimal.valueOf(k)), mathContext) : ZERO);
        M x = identidadInterna(filas);
        BigDecimal km1 = BigDecimal.valueOf(k - 1L), kk = BigDecimal.valueOf(k);
        for (int it = 0; it < 200; it++) {
            M siguiente = x.productoEscalar(km1).suma(producto(x.potencia(k - 1).inversa())).dividirPorEscalar(kk);
            if (siguiente.distanciaFrobenius(x).compareTo(tolerancia) <= 0) return siguiente;
            x = siguiente;
        }
        throw new comun.excepciones.ConvergenciaException("La raíz principal no convergió");
    }

    public final M seno() { return serieTrigonometrica(true, false); }
    public final M coseno() { return serieTrigonometrica(false, false); }
    public final M senoHiperbolico() { return serieTrigonometrica(true, true); }
    public final M cosenoHiperbolico() { return serieTrigonometrica(false, true); }
    private M serieTrigonometrica(boolean impar, boolean hiperbolica) {
        exigirCuadrada();
        M suma = impar ? comoEste(this) : identidadInterna(filas);
        M termino = comoEste(suma);
        int inicio = impar ? 1 : 0;
        M a2 = producto(this);
        for (int n = 1; n < 200; n++) {
            long d1 = inicio + 2L * n - 1, d2 = inicio + 2L * n;
            termino = termino.producto(a2).dividirPorEscalar(BigDecimal.valueOf(d1 * d2));
            if (!hiperbolica) termino = termino.opuesta();
            suma = suma.suma(termino);
            if (termino.normaFrobenius().compareTo(tolerancia) <= 0) return suma;
        }
        throw new comun.excepciones.ConvergenciaException("Serie matricial no convergente");
    }
    public final M tangente() { return seno().dividirPorDerecha(coseno()); }
    public final M signoMatricial() { return producto(potencia(2).raizCuadradaPrincipal().inversa()); }
    public final M valorAbsolutoMatricial() { return adjuntaHermitica().producto(this).raizCuadradaPrincipal(); }
    public final M funcionRacional(BigComplex[] numerador, BigComplex[] denominador) {
        return evaluarPolinomio(numerador).dividirPorDerecha(evaluarPolinomio(denominador));
    }
    public final M funcionAnalitica(UnaryOperator<BigComplex> funcion) { return funcionPorSchur(funcion); }
    public final M funcionPorJordan(UnaryOperator<BigComplex> funcion) {
        ResultadoDiagonalizacion<M> d = diagonalizacion();
        M fd = d.d().mapPosicional((i, j, z) -> i == j ? funcion.apply(z) : ZERO);
        return d.p().producto(fd).producto(d.inversaP());
    }
    public final M funcionPorSchur(UnaryOperator<BigComplex> funcion) {
        if (esNormal()) return funcionPorJordan(funcion);
        throw new OperacionNoDefinidaException("La evaluación escalar de Schur requiere matriz normal");
    }

    public final ResultadoIterativo<M> resolverJacobi(MatrizCompleja b, BigDecimal tol, int maxIteraciones) {
        exigirSistemaIterativo(b, tol, maxIteraciones);
        M bm = comoEste(b), x = crear(nuevosDatos(columnas, b.columnas(), ZERO), mathContext, tolerancia);
        BigDecimal residuo = residuoAXmenosB(x, bm);
        for (int it = 1; it <= maxIteraciones; it++) {
            BigComplex[][] siguiente = nuevosDatos(columnas, b.columnas(), ZERO);
            for (int i = 0; i < filas; i++) {
                if (cero(get(i, i))) throw new DivisionPorCeroException("Jacobi requiere diagonal no nula");
                for (int c = 0; c < b.columnas(); c++) {
                    BigComplex s = bm.get(i, c);
                    for (int j = 0; j < columnas; j++) if (j != i)
                        s = s.subtract(get(i, j).multiply(x.get(j, c), mathContext), mathContext);
                    siguiente[i][c] = s.divide(get(i, i), mathContext);
                }
            }
            x = crear(siguiente, mathContext, tolerancia);
            residuo = residuoAXmenosB(x, bm);
            if (residuo.compareTo(tol) <= 0) return new ResultadoIterativo<>(x, it, residuo, true);
        }
        return new ResultadoIterativo<>(x, maxIteraciones, residuo, false);
    }

    public final ResultadoIterativo<M> resolverGaussSeidel(MatrizCompleja b, BigDecimal tol, int maxIteraciones) {
        return resolverSOR(b, BigDecimal.ONE, tol, maxIteraciones);
    }

    public final ResultadoIterativo<M> resolverSOR(MatrizCompleja b, BigDecimal omega,
                                                   BigDecimal tol, int maxIteraciones) {
        exigirSistemaIterativo(b, tol, maxIteraciones);
        if (omega == null || omega.signum() <= 0 || omega.compareTo(BigDecimal.TWO) >= 0)
            throw new IllegalArgumentException("SOR exige 0 < omega < 2");
        M bm = comoEste(b);
        BigComplex[][] xd = nuevosDatos(columnas, b.columnas(), ZERO);
        BigDecimal residuo = new BigDecimal("1E100");
        for (int it = 1; it <= maxIteraciones; it++) {
            for (int i = 0; i < filas; i++) {
                if (cero(get(i, i))) throw new DivisionPorCeroException("SOR requiere diagonal no nula");
                for (int c = 0; c < b.columnas(); c++) {
                    BigComplex s = bm.get(i, c);
                    for (int j = 0; j < columnas; j++) if (j != i)
                        s = s.subtract(get(i, j).multiply(xd[j][c], mathContext), mathContext);
                    BigComplex gs = s.divide(get(i, i), mathContext);
                    xd[i][c] = xd[i][c].multiply(BigDecimal.ONE.subtract(omega), mathContext)
                            .add(gs.multiply(omega, mathContext), mathContext);
                }
            }
            M x = crear(xd, mathContext, tolerancia);
            residuo = residuoAXmenosB(x, bm);
            if (residuo.compareTo(tol) <= 0) return new ResultadoIterativo<>(x, it, residuo, true);
        }
        return new ResultadoIterativo<>(crear(xd, mathContext, tolerancia), maxIteraciones, residuo, false);
    }

    public final ResultadoIterativo<M> gradienteConjugado(MatrizCompleja b, BigDecimal tol, int maxIteraciones) {
        exigirSistemaIterativo(b, tol, maxIteraciones);
        if (!esHermiticaDefinidaPositiva()) throw new NoDefinidaPositivaException("CG exige matriz hermítica definida positiva");
        if (b.columnas() != 1) throw new DimensionMatrizException("CG implementa un vector derecho");
        M bm = comoEste(b), x = crear(nuevosDatos(columnas, 1, ZERO), mathContext, tolerancia);
        M rvec = bm.resta(producto(x)), p = comoEste(rvec);
        BigComplex rsold = rvec.productoFrobenius(rvec);
        BigDecimal residuo = rsold.abs(mathContext).sqrt(mathContext);
        for (int it = 1; it <= maxIteraciones; it++) {
            M ap = producto(p);
            BigComplex alpha = rsold.divide(p.productoFrobenius(ap), mathContext);
            x = x.suma(p.productoEscalar(alpha));
            rvec = rvec.resta(ap.productoEscalar(alpha));
            BigComplex rsnew = rvec.productoFrobenius(rvec);
            residuo = rsnew.abs(mathContext).sqrt(mathContext);
            if (residuo.compareTo(tol) <= 0) return new ResultadoIterativo<>(x, it, residuo, true);
            p = rvec.suma(p.productoEscalar(rsnew.divide(rsold, mathContext)));
            rsold = rsnew;
        }
        return new ResultadoIterativo<>(x, maxIteraciones, residuo, false);
    }

    /**
     * GMRES completo sin reinicio, con base de Arnoldi y mínimos cuadrados en
     * cada iteración.
     */
    public final ResultadoIterativo<M> GMRES(MatrizCompleja b, BigDecimal tol, int maxIteraciones) {
        exigirSistemaIterativo(b, tol, maxIteraciones);
        if (b.columnas() != 1) throw new DimensionMatrizException("GMRES implementa un vector derecho");
        M bm = comoEste(b), x0 = crear(nuevosDatos(columnas, 1, ZERO), mathContext, tolerancia);
        M r0 = bm.resta(producto(x0));
        BigDecimal beta = r0.normaFrobenius();
        if (beta.compareTo(tol) <= 0) return new ResultadoIterativo<>(x0, 0, beta, true);
        List<M> q = new ArrayList<>();
        q.add(r0.dividirPorEscalar(beta));
        BigComplex[][] h = nuevosDatos(maxIteraciones + 1, maxIteraciones, ZERO);
        M mejor = x0;
        BigDecimal residuo = beta;
        for (int k = 0; k < maxIteraciones; k++) {
            M v = producto(q.get(k));
            for (int j = 0; j <= k; j++) {
                h[j][k] = q.get(j).productoFrobenius(v);
                v = v.resta(q.get(j).productoEscalar(h[j][k]));
            }
            h[k + 1][k] = BigComplex.valueOf(v.normaFrobenius());
            if (!cero(h[k + 1][k])) q.add(v.dividirPorEscalar(h[k + 1][k]));
            BigComplex[][] hk = new BigComplex[k + 2][k + 1];
            for (int i = 0; i < k + 2; i++) System.arraycopy(h[i], 0, hk[i], 0, k + 1);
            M hm = crear(hk, mathContext, tolerancia);
            BigComplex[][] e = nuevosDatos(k + 2, 1, ZERO); e[0][0] = BigComplex.valueOf(beta);
            M y = hm.minimosCuadrados(crear(e, mathContext, tolerancia));
            M qmat = ensamblarColumnas(q.subList(0, k + 1));
            mejor = x0.suma(qmat.producto(y));
            residuo = residuoAXmenosB(mejor, bm);
            if (residuo.compareTo(tol) <= 0) return new ResultadoIterativo<>(mejor, k + 1, residuo, true);
            if (q.size() == k + 1) break;
        }
        return new ResultadoIterativo<>(mejor, Math.min(maxIteraciones, q.size()), residuo, false);
    }

    public final ResultadoIterativo<M> BiCG(MatrizCompleja b, BigDecimal tol, int maxIteraciones) {
        return BiCGSTAB(b, tol, maxIteraciones);
    }

    public final ResultadoIterativo<M> BiCGSTAB(MatrizCompleja b, BigDecimal tol, int maxIteraciones) {
        exigirSistemaIterativo(b, tol, maxIteraciones);
        if (b.columnas() != 1) throw new DimensionMatrizException("BiCGSTAB implementa un vector derecho");
        M bm = comoEste(b), x = crear(nuevosDatos(columnas, 1, ZERO), mathContext, tolerancia);
        M rv = bm.resta(producto(x)), r0 = comoEste(rv);
        M p = comoEste(rv), v = crear(nuevosDatos(filas, 1, ZERO), mathContext, tolerancia);
        BigComplex rhoOld = ONE, alpha = ONE, omega = ONE;
        BigDecimal residuo = rv.normaFrobenius();
        for (int it = 1; it <= maxIteraciones; it++) {
            BigComplex rho = r0.productoFrobenius(rv);
            if (cero(rho)) break;
            if (it > 1) {
                BigComplex beta = rho.divide(rhoOld, mathContext).multiply(alpha.divide(omega, mathContext), mathContext);
                p = rv.suma(p.resta(v.productoEscalar(omega)).productoEscalar(beta));
            }
            v = producto(p);
            alpha = rho.divide(r0.productoFrobenius(v), mathContext);
            M s = rv.resta(v.productoEscalar(alpha));
            if (s.normaFrobenius().compareTo(tol) <= 0) {
                x = x.suma(p.productoEscalar(alpha));
                return new ResultadoIterativo<>(x, it, s.normaFrobenius(), true);
            }
            M t = producto(s);
            omega = t.productoFrobenius(s).divide(t.productoFrobenius(t), mathContext);
            x = x.suma(p.productoEscalar(alpha)).suma(s.productoEscalar(omega));
            rv = s.resta(t.productoEscalar(omega));
            residuo = rv.normaFrobenius();
            if (residuo.compareTo(tol) <= 0) return new ResultadoIterativo<>(x, it, residuo, true);
            if (cero(omega)) break;
            rhoOld = rho;
        }
        return new ResultadoIterativo<>(x, maxIteraciones, residuo, false);
    }

    public final ResultadoAutovalor<M> metodoPotenciaInversa(BigDecimal tol, int maxIteraciones) {
        return inversa().metodoDeLaPotencia(tol, maxIteraciones);
    }
    public final ResultadoAutovalor<M> potenciaInversa(BigDecimal tol, int maxIteraciones) {
        ResultadoAutovalor<M> r = metodoPotenciaInversa(tol, maxIteraciones);
        return new ResultadoAutovalor<>(r.autovalor().reciprocal(mathContext), r.autovector(),
                r.iteraciones(), r.residuo(), r.convergio());
    }
    public final ResultadoAutovalor<M> iteracionCocienteRayleigh(BigDecimal tol, int maxIteraciones) {
        exigirCuadrada();
        M v = crear(nuevosDatos(filas, 1, ONE), mathContext, tolerancia).normalizarColumnas();
        BigComplex mu = v.adjuntaHermitica().producto(producto(v)).get(0, 0);
        BigDecimal residuo = new BigDecimal("1E100");
        for (int it = 1; it <= maxIteraciones; it++) {
            M desplazada = resta(identidadInterna(filas).productoEscalar(mu));
            M y;
            try { y = desplazada.resolver(v); }
            catch (MatrizSingularException | OperacionNoDefinidaException e) {
                return new ResultadoAutovalor<>(mu, v, it, BigDecimal.ZERO, true);
            }
            v = y.normalizarColumnas();
            mu = v.adjuntaHermitica().producto(producto(v)).get(0, 0);
            residuo = producto(v).resta(v.productoEscalar(mu)).normaFrobenius();
            if (residuo.compareTo(tol) <= 0) return new ResultadoAutovalor<>(mu, v, it, residuo, true);
        }
        return new ResultadoAutovalor<>(mu, v, maxIteraciones, residuo, false);
    }
    public final ResultadoQR<M> algoritmoQR(BigDecimal tol, int maxIteraciones) {
        M t = algoritmoQRParaAutovalores(tol, maxIteraciones);
        return new ResultadoQR<>(identidadInterna(filas), t, t.resta(t.triangularSuperior()).normaFrobenius(), BigDecimal.ZERO);
    }
    public final ResultadoQR<M> Lanczos(int pasos) {
        if (!esHermitica()) throw new MatrizNoHermiticaException("Lanczos exige hermiticidad");
        return Arnoldi(pasos);
    }
    public final ResultadoQR<M> Arnoldi(int pasos) {
        exigirCuadrada();
        if (pasos <= 0 || pasos > filas) throw new DimensionMatrizException("Pasos de Krylov inválidos");
        BigComplex[][] inicial = nuevosDatos(filas, 1, ONE);
        List<M> q = new ArrayList<>();
        q.add(crear(inicial, mathContext, tolerancia).normalizarColumnas());
        BigComplex[][] h = nuevosDatos(pasos, pasos, ZERO);
        for (int k = 0; k < pasos; k++) {
            M v = producto(q.get(k));
            for (int j = 0; j <= k; j++) {
                h[j][k] = q.get(j).productoFrobenius(v);
                v = v.resta(q.get(j).productoEscalar(h[j][k]));
            }
            if (k + 1 < pasos) {
                h[k + 1][k] = BigComplex.valueOf(v.normaFrobenius());
                if (cero(h[k + 1][k])) break;
                q.add(v.dividirPorEscalar(h[k + 1][k]));
            }
        }
        M qm = ensamblarColumnas(q), hm = crear(h, mathContext, tolerancia);
        return new ResultadoQR<>(qm, hm, producto(qm).resta(qm.producto(hm)).normaFrobenius(),
                qm.adjuntaHermitica().producto(qm).distanciaFrobenius(identidadInterna(qm.columnas())));
    }
    public final ResultadoIterativo<M> refinamientoIterativo(MatrizCompleja b, BigDecimal tol, int maxIteraciones) {
        M bm = comoEste(b), x = resolver(b);
        BigDecimal residuo = residuoAXmenosB(x, bm);
        for (int it = 1; it <= maxIteraciones; it++) {
            if (residuo.compareTo(tol) <= 0) return new ResultadoIterativo<>(x, it - 1, residuo, true);
            M correccion = resolver(bm.resta(producto(x)));
            x = x.suma(correccion);
            residuo = residuoAXmenosB(x, bm);
        }
        return new ResultadoIterativo<>(x, maxIteraciones, residuo, false);
    }
    public final M precondicionarJacobi() {
        exigirCuadrada();
        BigComplex[][] d = nuevosDatos(filas, filas, ZERO);
        for (int i = 0; i < filas; i++) {
            if (cero(get(i, i))) throw new DivisionPorCeroException("Diagonal nula");
            d[i][i] = get(i, i).reciprocal(mathContext);
        }
        return crear(d, mathContext, tolerancia);
    }
    public final M precondicionarLUIncompleta() { return descomposicionLUP().l(); }

    private void exigirSistemaIterativo(MatrizCompleja b, BigDecimal tol, int max) {
        exigirCuadrada(); mismaPolitica(b);
        if (b.filas() != filas) throw new DimensionMatrizException("Término independiente incompatible");
        Complejos.validarTolerancia(tol);
        if (max <= 0) throw new IllegalArgumentException("Máximo de iteraciones no positivo");
    }

    public final List<List<M>> particionar(int[] cortesFilas, int[] cortesColumnas) {
        validarParticion(cortesFilas, filas); validarParticion(cortesColumnas, columnas);
        int[] rf = limites(cortesFilas, filas), rc = limites(cortesColumnas, columnas);
        List<List<M>> salida = new ArrayList<>();
        for (int i = 0; i + 1 < rf.length; i++) {
            List<M> filaBloques = new ArrayList<>();
            for (int j = 0; j + 1 < rc.length; j++)
                filaBloques.add(submatriz(rf[i], rf[i + 1], rc[j], rc[j + 1]));
            salida.add(List.copyOf(filaBloques));
        }
        return List.copyOf(salida);
    }

    public final M obtenerBloque(int filaDesde, int filaHasta, int columnaDesde, int columnaHasta) {
        return submatriz(filaDesde, filaHasta, columnaDesde, columnaHasta);
    }

    public final M ensamblarBloques(MatrizCompleja[][] bloques) {
        Objects.requireNonNull(bloques);
        if (bloques.length == 0 || bloques[0] == null || bloques[0].length == 0)
            throw new DimensionMatrizException("Partición vacía");
        int bc = bloques[0].length;
        int totalF = 0, totalC = 0;
        int[] altos = new int[bloques.length], anchos = new int[bc];
        for (int i = 0; i < bloques.length; i++) {
            if (bloques[i] == null || bloques[i].length != bc) throw new DimensionMatrizException("Bloques irregulares");
            altos[i] = Objects.requireNonNull(bloques[i][0]).filas();
            totalF += altos[i];
            for (int j = 0; j < bc; j++) {
                MatrizCompleja b = Objects.requireNonNull(bloques[i][j]);
                mismaPolitica(b);
                if (b.filas() != altos[i]) throw new DimensionMatrizException("Alturas de bloque incompatibles");
                if (i == 0) { anchos[j] = b.columnas(); totalC += anchos[j]; }
                else if (b.columnas() != anchos[j]) throw new DimensionMatrizException("Anchuras de bloque incompatibles");
            }
        }
        BigComplex[][] d = new BigComplex[totalF][totalC];
        int fi = 0;
        for (int bi = 0; bi < bloques.length; bi++) {
            int cj = 0;
            for (int bj = 0; bj < bc; bj++) {
                MatrizCompleja b = bloques[bi][bj];
                for (int i = 0; i < b.filas(); i++) for (int j = 0; j < b.columnas(); j++) d[fi + i][cj + j] = b.get(i, j);
                cj += b.columnas();
            }
            fi += altos[bi];
        }
        return crear(d, mathContext, tolerancia);
    }

    public final M sumaPorBloques(MatrizCompleja otra, int[] cortesFilas, int[] cortesColumnas) {
        esMatrizPorBloques(cortesFilas, cortesColumnas);
        return suma(otra);
    }
    public final M productoPorBloques(MatrizCompleja otra, int corteComun) {
        if (corteComun <= 0 || corteComun >= columnas || columnas != otra.filas())
            throw new DimensionMatrizException("Corte de producto por bloques inválido");
        return producto(otra);
    }
    public final M traspuestaPorBloques(int[] cortesFilas, int[] cortesColumnas) {
        esMatrizPorBloques(cortesFilas, cortesColumnas); return traspuesta();
    }
    public final M adjuntaPorBloques(int[] cortesFilas, int[] cortesColumnas) {
        esMatrizPorBloques(cortesFilas, cortesColumnas); return adjuntaHermitica();
    }
    public final M eliminacionGaussPorBloques(int corte) {
        complementoDeSchur(corte); return formaEscalonada();
    }
    public final M diagonalPorBloques(List<? extends MatrizCompleja> bloques) {
        if (bloques == null || bloques.isEmpty()) throw new DimensionMatrizException("Sin bloques");
        int totalF = bloques.stream().mapToInt(MatrizCompleja::filas).sum();
        int totalC = bloques.stream().mapToInt(MatrizCompleja::columnas).sum();
        BigComplex[][] d = nuevosDatos(totalF, totalC, ZERO);
        int f = 0, c = 0;
        for (MatrizCompleja b : bloques) {
            mismaPolitica(b);
            for (int i = 0; i < b.filas(); i++) for (int j = 0; j < b.columnas(); j++) d[f + i][c + j] = b.get(i, j);
            f += b.filas(); c += b.columnas();
        }
        return crear(d, mathContext, tolerancia);
    }
    public final M resolverSistemaPorBloques(MatrizCompleja b, int corte) {
        complementoDeSchur(corte); return resolver(b);
    }

    private static void validarParticion(int[] cortes, int limite) {
        Objects.requireNonNull(cortes);
        int previo = 0;
        for (int c : cortes) {
            if (c <= previo || c >= limite) throw new DimensionMatrizException("Cortes no crecientes o fuera de rango");
            previo = c;
        }
    }
    private static void validarTamanosBloques(int[] tamanos, int total) {
        if (tamanos == null || tamanos.length == 0 || Arrays.stream(tamanos).anyMatch(x -> x <= 0)
                || Arrays.stream(tamanos).sum() != total) throw new DimensionMatrizException("Tamaños de bloque inválidos");
    }
    private static int[] limites(int[] cortes, int limite) {
        int[] r = new int[cortes.length + 2];
        r[0] = 0; System.arraycopy(cortes, 0, r, 1, cortes.length); r[r.length - 1] = limite;
        return r;
    }

    public final M centrarColumnas() {
        BigComplex[][] d = copiaDatos();
        for (int j = 0; j < columnas; j++) {
            BigComplex media = ZERO;
            for (int i = 0; i < filas; i++) media = media.add(get(i, j), mathContext);
            media = media.divide(BigDecimal.valueOf(filas), mathContext);
            for (int i = 0; i < filas; i++) d[i][j] = get(i, j).subtract(media, mathContext);
        }
        return crear(d, mathContext, tolerancia);
    }
    public final M centrarFilas() { return traspuesta().centrarColumnas().traspuesta(); }
    public final M estandarizarColumnas() {
        M c = centrarColumnas();
        BigComplex[][] d = c.copiaDatos();
        for (int j = 0; j < columnas; j++) {
            BigDecimal s = BigDecimal.ZERO;
            for (int i = 0; i < filas; i++) s = s.add(c.get(i, j).absSquare(mathContext), mathContext);
            BigDecimal desv = s.divide(BigDecimal.valueOf(Math.max(1, filas - 1)), mathContext).sqrt(mathContext);
            if (desv.compareTo(tolerancia) <= 0) throw new DivisionPorCeroException("Columna constante");
            for (int i = 0; i < filas; i++) d[i][j] = d[i][j].divide(desv, mathContext);
        }
        return crear(d, mathContext, tolerancia);
    }
    public final M matrizCovarianza() {
        if (filas < 2) throw new DimensionMatrizException("Covarianza requiere dos observaciones");
        M c = centrarColumnas();
        return c.adjuntaHermitica().producto(c).dividirPorEscalar(BigDecimal.valueOf(filas - 1L));
    }
    public final M matrizDispersion() { M c = centrarColumnas(); return c.adjuntaHermitica().producto(c); }
    public final M matrizCorrelacion() {
        M z = estandarizarColumnas();
        return z.adjuntaHermitica().producto(z).dividirPorEscalar(BigDecimal.valueOf(filas - 1L));
    }
    public final ResultadoDiagonalizacion<M> PCA() { return matrizCovarianza().diagonalizacionUnitaria(); }
    public final M blanqueamiento() {
        ResultadoDiagonalizacion<M> pca = PCA();
        M invSqrt = pca.d().mapPosicional((i, j, z) -> i == j
                ? (cero(z) ? ZERO : BigComplexMath.sqrt(z, mathContext).reciprocal(mathContext)) : ZERO);
        return centrarColumnas().producto(pca.p()).producto(invSqrt);
    }
    public final M regresionMinimosCuadrados(MatrizCompleja y) { return minimosCuadrados(y); }
    public final M matrizSombrero() { return producto(pseudoinversaMoorePenrose()); }

    public final M matrizConfusion(int[] reales, int[] predichos, int clases) {
        if (reales == null || predichos == null || reales.length != predichos.length || clases <= 0)
            throw new DimensionMatrizException("Etiquetas inválidas");
        BigComplex[][] d = nuevosDatos(clases, clases, ZERO);
        for (int i = 0; i < reales.length; i++) {
            if (reales[i] < 0 || reales[i] >= clases || predichos[i] < 0 || predichos[i] >= clases)
                throw new IndiceMatrizException("Clase fuera de rango");
            d[reales[i]][predichos[i]] = d[reales[i]][predichos[i]].add(ONE);
        }
        return crear(d, mathContext, tolerancia);
    }

    public final M matrizGrados() {
        exigirAdyacencia();
        BigComplex[][] d = nuevosDatos(filas, filas, ZERO);
        for (int i = 0; i < filas; i++) {
            BigComplex s = ZERO;
            for (int j = 0; j < columnas; j++) s = s.add(get(i, j), mathContext);
            d[i][i] = s;
        }
        return crear(d, mathContext, tolerancia);
    }
    public final M laplaciana() { return matrizGrados().resta(this); }
    public final M laplacianaNormalizada() {
        M d = matrizGrados();
        BigComplex[][] inv = nuevosDatos(filas, filas, ZERO);
        for (int i = 0; i < filas; i++)
            if (!cero(d.get(i, i))) inv[i][i] = BigComplexMath.sqrt(d.get(i, i), mathContext).reciprocal(mathContext);
        M di = crear(inv, mathContext, tolerancia);
        return identidadInterna(filas).resta(di.producto(this).producto(di));
    }
    public final M numeroDeRecorridos(int longitud) {
        if (longitud < 0) throw new IllegalArgumentException("Longitud negativa");
        exigirAdyacencia(); return potencia(longitud);
    }
    public final boolean[][] clausuraTransitivaBooleana() {
        exigirAdyacencia();
        boolean[][] r = new boolean[filas][filas];
        for (int i = 0; i < filas; i++) for (int j = 0; j < filas; j++) r[i][j] = i == j || !cero(get(i, j));
        for (int k = 0; k < filas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < filas; j++)
            r[i][j] |= r[i][k] && r[k][j];
        return r;
    }
    public final boolean esGrafoConexo() { return componentesConexas().size() == 1; }
    public final List<List<Integer>> componentesConexas() {
        exigirCuadrada();
        boolean[] visto = new boolean[filas];
        List<List<Integer>> comps = new ArrayList<>();
        for (int s = 0; s < filas; s++) if (!visto[s]) {
            List<Integer> comp = new ArrayList<>(), cola = new ArrayList<>(); cola.add(s); visto[s] = true;
            for (int p = 0; p < cola.size(); p++) {
                int i = cola.get(p); comp.add(i);
                for (int j = 0; j < filas; j++) if (!visto[j] && (!cero(get(i, j)) || !cero(get(j, i)))) {
                    visto[j] = true; cola.add(j);
                }
            }
            comps.add(List.copyOf(comp));
        }
        return List.copyOf(comps);
    }
    private void exigirAdyacencia() {
        exigirCuadrada();
        if (!esNoNegativa()) throw new OperacionNoDefinidaException("La adyacencia debe ser real no negativa");
    }

    public final M distribuirTrasPasos(BigComplex[] distribucionFila, int pasos) {
        if (!esMatrizDeTransicion()) throw new OperacionNoDefinidaException("Se requiere matriz de transición");
        BigComplex[] r = potencia(pasos).vectorFilaPorMatriz(distribucionFila);
        return crear(new BigComplex[][]{r}, mathContext, tolerancia);
    }
    public final M distribucionEstacionaria() {
        if (!esMatrizDeTransicion()) throw new OperacionNoDefinidaException("Se requiere matriz de transición");
        M sistema = traspuesta().resta(identidadInterna(filas));
        BigComplex[][] filaUnos = new BigComplex[1][filas]; Arrays.fill(filaUnos[0], ONE);
        M a = sistema.eliminarFila(filas - 1).concatenarVertical(crear(filaUnos, mathContext, tolerancia));
        BigComplex[][] b = nuevosDatos(filas, 1, ZERO); b[filas - 1][0] = ONE;
        return a.resolver(crear(b, mathContext, tolerancia)).traspuesta();
    }
    public final boolean esCadenaRegular() {
        if (!esMatrizDeTransicion()) return false;
        M p = comoEste(this);
        for (int k = 1; k <= filas * filas; k++) {
            if (p.esPositivaPorEntradas()) return true;
            p = p.producto(this);
        }
        return false;
    }

    public final M vectorizarPorColumnas() {
        BigComplex[][] d = new BigComplex[filas * columnas][1];
        int k = 0;
        for (int j = 0; j < columnas; j++) for (int i = 0; i < filas; i++) d[k++][0] = get(i, j);
        return crear(d, mathContext, tolerancia);
    }
    public final M vectorizarPorFilas() {
        BigComplex[][] d = new BigComplex[filas * columnas][1];
        int k = 0;
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++) d[k++][0] = get(i, j);
        return crear(d, mathContext, tolerancia);
    }
    public final M desvectorizar(int nuevasFilas, int nuevasColumnas, boolean porColumnas) {
        if (columnas != 1 || nuevasFilas * nuevasColumnas != filas)
            throw new DimensionMatrizException("Vector y forma incompatibles");
        BigComplex[][] d = new BigComplex[nuevasFilas][nuevasColumnas];
        int k = 0;
        if (porColumnas) for (int j = 0; j < nuevasColumnas; j++) for (int i = 0; i < nuevasFilas; i++) d[i][j] = get(k++, 0);
        else for (int i = 0; i < nuevasFilas; i++) for (int j = 0; j < nuevasColumnas; j++) d[i][j] = get(k++, 0);
        return crear(d, mathContext, tolerancia);
    }
    public final M sumaKronecker(MatrizCompleja b) {
        exigirCuadrada();
        M bm = comoEste(b); bm.exigirCuadrada();
        return productoKronecker(bm.identidadInterna(bm.filas())).suma(identidadInterna(filas).productoKronecker(bm));
    }
    public final boolean verificarIdentidadVecAXB(MatrizCompleja x, MatrizCompleja b) {
        M xm = comoEste(x), bm = comoEste(b);
        M izquierda = producto(xm).producto(bm).vectorizarPorColumnas();
        M derecha = bm.traspuesta().productoKronecker(this).producto(xm.vectorizarPorColumnas());
        return izquierda.aproximadamenteIgual(derecha, tolerancia);
    }
    public final M reshape(int nuevasFilas, int nuevasColumnas) {
        if (nuevasFilas * nuevasColumnas != filas * columnas) throw new DimensionMatrizException("Número de elementos distinto");
        return vectorizarPorFilas().desvectorizar(nuevasFilas, nuevasColumnas, false);
    }
    public final BigComplex[] flatten() {
        BigComplex[] salida = new BigComplex[filas * columnas];
        int k = 0;
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++) salida[k++] = get(i, j);
        return salida;
    }
    public final M permutarIndices(int[] permutacionFilas, int[] permutacionColumnas) {
        validarPermutacion(permutacionFilas, filas); validarPermutacion(permutacionColumnas, columnas);
        return extraerFilas(permutacionFilas).extraerColumnas(permutacionColumnas);
    }
    private static void validarPermutacion(int[] p, int n) {
        if (p == null || p.length != n || Arrays.stream(p).distinct().count() != n
                || Arrays.stream(p).anyMatch(x -> x < 0 || x >= n)) throw new DimensionMatrizException("Permutación inválida");
    }

    public final boolean verificarInversa(MatrizCompleja candidata, BigDecimal tol) {
        M c = comoEste(candidata);
        return esCuadrada() && c.filas() == filas && producto(c).aproximadamenteIgual(identidadInterna(filas), tol)
                && c.producto(this).aproximadamenteIgual(identidadInterna(filas), tol);
    }
    public final boolean verificarLU(ResultadoLU<M> lu, BigDecimal tol) {
        return lu.p().producto(this).aproximadamenteIgual(lu.l().producto(lu.u()), tol);
    }
    public final boolean verificarQR(ResultadoQR<M> qr, BigDecimal tol) {
        return aproximadamenteIgual(qr.q().producto(qr.r()), tol)
                && qr.q().adjuntaHermitica().producto(qr.q()).aproximadamenteIgual(identidadInterna(qr.q().columnas()), tol);
    }
    public final boolean verificarCholesky(ResultadoLU<M> ch, BigDecimal tol) {
        return aproximadamenteIgual(ch.l().producto(ch.l().adjuntaHermitica()), tol);
    }
    public final boolean verificarSVD(ResultadoSVD<M> svd, BigDecimal tol) {
        return aproximadamenteIgual(svd.u().producto(svd.sigma()).producto(svd.v().adjuntaHermitica()), tol);
    }
    public final boolean verificarDiagonalizacion(ResultadoDiagonalizacion<M> d, BigDecimal tol) {
        return aproximadamenteIgual(d.p().producto(d.d()).producto(d.inversaP()), tol);
    }
    public final boolean verificarSchur(ResultadoQR<M> schur, BigDecimal tol) {
        return aproximadamenteIgual(schur.q().producto(schur.r()).producto(schur.q().adjuntaHermitica()), tol);
    }
    public final boolean verificarFormaJordan(ResultadoDiagonalizacion<M> jordan, BigDecimal tol) {
        return verificarDiagonalizacion(jordan, tol);
    }
    public final boolean verificarPseudoinversa(MatrizCompleja p, BigDecimal tol) {
        M pm = comoEste(p);
        return producto(pm).producto(this).aproximadamenteIgual(this, tol)
                && pm.producto(this).producto(pm).aproximadamenteIgual(pm, tol)
                && producto(pm).esHermitica()
                && pm.producto(this).esHermitica();
    }
    public final boolean verificarSolucionSistema(MatrizCompleja x, MatrizCompleja b, BigDecimal tol) {
        return producto(x).aproximadamenteIgual(b, tol);
    }
    public final boolean verificarAutopar(BigComplex lambda, MatrizCompleja v, BigDecimal tol) {
        M vm = comoEste(v);
        return producto(vm).aproximadamenteIgual(vm.productoEscalar(lambda), tol);
    }
    public final boolean verificarProyector(BigDecimal tol) { return potencia(2).aproximadamenteIgual(this, tol); }
    public final boolean verificarUnitaria(BigDecimal tol) { return esUnitariaAproximadamente(tol); }
    public final boolean verificarHermitica(BigDecimal tol) { return aproximadamenteIgual(adjuntaHermitica(), tol); }
    public final boolean verificarCayleyHamilton(BigDecimal tol) {
        return evaluarPolinomio(coeficientesPolinomioCaracteristico()).normaFrobenius().compareTo(tol) <= 0;
    }
    public final BigDecimal residuoAXmenosB(MatrizCompleja x, MatrizCompleja b) {
        return producto(x).resta(b).normaFrobenius();
    }
    public final BigDecimal errorReconstruccion(MatrizCompleja reconstruida) {
        return distanciaFrobenius(reconstruida);
    }

    public final M derivarElementoAElemento(UnaryOperator<BigComplex> derivada) { return map(derivada); }
    public final M integrarElementoAElemento(UnaryOperator<BigComplex> primitiva) { return map(primitiva); }
    public final M diferencialDeProducto(MatrizCompleja dA, MatrizCompleja b, MatrizCompleja dB) {
        M da = comoEste(dA), bm = comoEste(b), db = comoEste(dB);
        return da.producto(bm).suma(producto(db));
    }
    public final M diferencialDeInversa(MatrizCompleja dA) {
        M inv = inversa();
        return inv.producto(dA).producto(inv).opuesta();
    }
    public final BigComplex diferencialDeDeterminante(MatrizCompleja dA) {
        return determinante().multiply(inversa().producto(dA).traza(), mathContext);
    }
    public final M gradienteDeTraza() { exigirCuadrada(); return identidadInterna(filas); }
    public final M gradienteFormaCuadratica(BigComplex[] x) {
        if (x.length != filas) throw new DimensionMatrizException("Vector incompatible");
        BigComplex[][] d = new BigComplex[filas][filas];
        for (int i = 0; i < filas; i++) for (int j = 0; j < filas; j++)
            d[i][j] = x[i].conjugate().multiply(x[j], mathContext);
        return crear(d, mathContext, tolerancia);
    }
    public final M jacobiano(UnaryOperator<BigComplex> derivada) { return derivarElementoAElemento(derivada); }
    public final M hessiano(UnaryOperator<BigComplex> segundaDerivada) { return map(segundaDerivada); }

    public final M derivadaNumerica(java.util.function.Function<BigDecimal, ? extends MatrizCompleja> funcion,
                                    BigDecimal punto, BigDecimal paso) {
        Objects.requireNonNull(funcion); Objects.requireNonNull(punto); Objects.requireNonNull(paso);
        if (paso.signum() == 0) throw new DivisionPorCeroException("Paso nulo");
        M mas = comoEste(funcion.apply(punto.add(paso, mathContext)));
        M menos = comoEste(funcion.apply(punto.subtract(paso, mathContext)));
        return mas.resta(menos).dividirPorEscalar(paso.multiply(BigDecimal.TWO, mathContext));
    }
    public final M integralNumerica(java.util.function.Function<BigDecimal, ? extends MatrizCompleja> funcion,
                                   BigDecimal a, BigDecimal b, int intervalos) {
        Objects.requireNonNull(funcion);
        if (intervalos <= 0) throw new IllegalArgumentException("Intervalos no positivos");
        BigDecimal h = b.subtract(a, mathContext).divide(BigDecimal.valueOf(intervalos), mathContext);
        M suma = comoEste(funcion.apply(a)).suma(comoEste(funcion.apply(b))).dividirPorEscalar(BigDecimal.TWO);
        for (int i = 1; i < intervalos; i++) suma = suma.suma(comoEste(funcion.apply(a.add(h.multiply(BigDecimal.valueOf(i)), mathContext))));
        return suma.productoEscalar(h);
    }

    public final String formato(int escala) {
        if (escala < 0) throw new IllegalArgumentException("Escala negativa");
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < filas; i++) {
            sb.append('[');
            for (int j = 0; j < columnas; j++) {
                if (j > 0) sb.append(", ");
                BigComplex z = get(i, j);
                sb.append(z.re.setScale(escala, RoundingMode.HALF_EVEN));
                if (z.im.signum() >= 0) sb.append('+');
                sb.append(z.im.setScale(escala, RoundingMode.HALF_EVEN)).append('i');
            }
            sb.append(']');
            if (i + 1 < filas) sb.append(System.lineSeparator());
        }
        return sb.toString();
    }
    public final String serializar() {
        StringBuilder sb = new StringBuilder(filas + "x" + columnas + ":");
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++) {
            if (i != 0 || j != 0) sb.append(';');
            sb.append(get(i, j).re.toPlainString()).append(',').append(get(i, j).im.toPlainString());
        }
        return sb.toString();
    }

    @Override public final boolean equals(Object objeto) {
        return objeto instanceof MatrizCompleja otra && igualExacta(otra);
    }
    @Override public final int hashCode() {
        int h = Objects.hash(filas, columnas);
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++) h = 31 * h + get(i, j).hashCode();
        return h;
    }
    @Override public final String toString() { return formato(Math.min(8, mathContext.getPrecision())); }
}
