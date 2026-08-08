package comun;

import ch.obermuhlner.math.big.BigComplex;
import ch.obermuhlner.math.big.BigDecimalMath;
import comun.excepciones.DimensionMatrizException;
import comun.excepciones.DivisionPorCeroException;
import comun.excepciones.IndiceMatrizException;
import comun.excepciones.MatrizException;
import comun.excepciones.OperacionNoDefinidaException;
import comun.resultados.ResultadoHOSVD;
import comun.resultados.ResultadoDiagonalizacion;
import comun.resultados.ResultadoSVD;
import comun.resultados.ResultadoTSVD;

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

    /** Devuelve el plano perpendicular a {@code eje} en {@code indice}. */
    public final P plano(int eje, int indice) {
        validarEje(eje);
        return switch (eje) {
            case 0 -> capa(indice);
            case 1 -> cortePorFila(indice);
            default -> cortePorColumna(indice);
        };
    }

    /**
     * Extrae una fibra variando un eje y fijando los otros dos. Para el eje 0
     * los índices son fila/columna; para el 1, capa/columna; para el 2,
     * capa/fila.
     */
    public final BigComplex[] fibra(int eje, int indiceA, int indiceB) {
        validarEje(eje);
        int longitud = tamanoEje(eje);
        BigComplex[] salida = new BigComplex[longitud];
        for (int n = 0; n < longitud; n++) {
            int[] p = coordenadasFibra(eje, n, indiceA, indiceB);
            validarIndice(p[0], p[1], p[2]);
            salida[n] = get(p[0], p[1], p[2]);
        }
        return salida;
    }

    public final T conFibra(int eje, int indiceA, int indiceB, BigComplex[] valores) {
        validarEje(eje);
        Objects.requireNonNull(valores, "La fibra no puede ser null");
        if (valores.length != tamanoEje(eje)) throw new DimensionMatrizException("Longitud de fibra incorrecta");
        BigComplex[][][] datos = copiaDatos();
        for (int n = 0; n < valores.length; n++) {
            int[] p = coordenadasFibra(eje, n, indiceA, indiceB);
            validarIndice(p[0], p[1], p[2]);
            datos[p[0]][p[1]][p[2]] = redondear(valores[n]);
        }
        return crear(datos, mathContext, tolerancia);
    }

    public final int tamanoEje(int eje) {
        validarEje(eje);
        return switch (eje) { case 0 -> capas; case 1 -> filas; default -> columnas; };
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

    /** Extrae índices, permitiendo repetición y reordenación, sobre cualquier eje. */
    public final T extraerIndices(int eje, int... indices) {
        validarEje(eje);
        Objects.requireNonNull(indices, "Los índices no pueden ser null");
        if (indices.length == 0) throw new DimensionMatrizException("Debe extraerse algún índice");
        int[] forma = {capas, filas, columnas};
        for (int indice : indices) if (indice < 0 || indice >= forma[eje])
            throw new IndiceMatrizException("Índice fuera del eje " + eje + ": " + indice);
        int[] salida = forma.clone();
        salida[eje] = indices.length;
        BigComplex[][][] datos = new BigComplex[salida[0]][salida[1]][salida[2]];
        for (int k = 0; k < salida[0]; k++) for (int i = 0; i < salida[1]; i++) for (int j = 0; j < salida[2]; j++) {
            int[] p = {k, i, j};
            p[eje] = indices[p[eje]];
            datos[k][i][j] = get(p[0], p[1], p[2]);
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

    public final T eliminarIndice(int eje, int indice) {
        validarEje(eje);
        int longitud = tamanoEje(eje);
        if (longitud == 1) throw new DimensionMatrizException("No se admite un eje vacío");
        if (indice < 0 || indice >= longitud) throw new IndiceMatrizException("Índice fuera del eje");
        int[] indices = new int[longitud - 1];
        for (int n = 0, p = 0; n < longitud; n++) if (n != indice) indices[p++] = n;
        return extraerIndices(eje, indices);
    }

    /** Inserta un plano perpendicular al eje indicado. */
    public final T insertarPlano(int eje, int posicion, MatrizCompleja plano) {
        validarEje(eje);
        Objects.requireNonNull(plano, "El plano no puede ser null");
        if (!mathContext.equals(plano.mathContext()) || tolerancia.compareTo(plano.tolerancia()) != 0)
            throw new OperacionNoDefinidaException("No se mezclan políticas de precisión distintas");
        int[] forma = {capas, filas, columnas};
        if (posicion < 0 || posicion > forma[eje]) throw new IndiceMatrizException("Posición fuera del eje");
        int esperadoFilas = forma[(eje + 1) % 3];
        int esperadoColumnas = forma[(eje + 2) % 3];
        if (eje == 1) { esperadoFilas = capas; esperadoColumnas = columnas; }
        if (plano.filas() != esperadoFilas || plano.columnas() != esperadoColumnas)
            throw new DimensionMatrizException("Dimensiones del plano incompatibles");
        int[] salida = forma.clone(); salida[eje]++;
        BigComplex[][][] datos = new BigComplex[salida[0]][salida[1]][salida[2]];
        for (int k = 0; k < salida[0]; k++) for (int i = 0; i < salida[1]; i++) for (int j = 0; j < salida[2]; j++) {
            int[] p = {k, i, j};
            if (p[eje] == posicion) {
                int[] q = coordenadasPlano(eje, k, i, j);
                datos[k][i][j] = plano.get(q[0], q[1]);
            } else {
                if (p[eje] > posicion) p[eje]--;
                datos[k][i][j] = get(p[0], p[1], p[2]);
            }
        }
        return crear(datos, mathContext, tolerancia);
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

    /** Divide por cortes internos estrictamente crecientes sobre un eje. */
    public final List<T> particionar(int eje, int... cortes) {
        validarEje(eje);
        Objects.requireNonNull(cortes, "Los cortes no pueden ser null");
        int anterior = 0;
        List<T> salida = new ArrayList<>(cortes.length + 1);
        for (int corte : cortes) {
            if (corte <= anterior || corte >= tamanoEje(eje))
                throw new IndiceMatrizException("Cortes no crecientes o fuera del eje");
            salida.add(intervaloEje(eje, anterior, corte));
            anterior = corte;
        }
        salida.add(intervaloEje(eje, anterior, tamanoEje(eje)));
        return List.copyOf(salida);
    }

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

    /** Expande únicamente ejes de longitud uno. */
    public final T broadcast(int nuevasCapas, int nuevasFilas, int nuevasColumnas) {
        int[] origen = {capas, filas, columnas};
        int[] salida = {nuevasCapas, nuevasFilas, nuevasColumnas};
        for (int eje = 0; eje < 3; eje++)
            if (salida[eje] <= 0 || (origen[eje] != salida[eje] && origen[eje] != 1))
                throw new DimensionMatrizException("Broadcast incompatible en el eje " + eje);
        BigComplex[][][] datos = new BigComplex[nuevasCapas][nuevasFilas][nuevasColumnas];
        for (int k = 0; k < nuevasCapas; k++) for (int i = 0; i < nuevasFilas; i++) for (int j = 0; j < nuevasColumnas; j++)
            datos[k][i][j] = get(capas == 1 ? 0 : k, filas == 1 ? 0 : i, columnas == 1 ? 0 : j);
        return crear(datos, mathContext, tolerancia);
    }

    /** Repite el tensor completo el número indicado de veces en cada eje. */
    public final T repetir(int repeticionesCapas, int repeticionesFilas, int repeticionesColumnas) {
        if (repeticionesCapas <= 0 || repeticionesFilas <= 0 || repeticionesColumnas <= 0)
            throw new DimensionMatrizException("Las repeticiones deben ser positivas");
        int nc = Math.multiplyExact(capas, repeticionesCapas);
        int nf = Math.multiplyExact(filas, repeticionesFilas);
        int nn = Math.multiplyExact(columnas, repeticionesColumnas);
        BigComplex[][][] datos = new BigComplex[nc][nf][nn];
        for (int k = 0; k < nc; k++) for (int i = 0; i < nf; i++) for (int j = 0; j < nn; j++)
            datos[k][i][j] = get(k % capas, i % filas, j % columnas);
        return crear(datos, mathContext, tolerancia);
    }

    /** Añade bordes constantes antes y después de cada eje. */
    public final T rellenarBordes(int capasAntes, int capasDespues, int filasAntes, int filasDespues,
                                  int columnasAntes, int columnasDespues, BigComplex valor) {
        int[] bordes = {capasAntes, capasDespues, filasAntes, filasDespues, columnasAntes, columnasDespues};
        if (Arrays.stream(bordes).anyMatch(x -> x < 0)) throw new DimensionMatrizException("Borde negativo");
        BigComplex[][][] datos = nuevosDatos(Math.addExact(capas, capasAntes + capasDespues),
                Math.addExact(filas, filasAntes + filasDespues),
                Math.addExact(columnas, columnasAntes + columnasDespues), redondear(valor));
        for (int k = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            datos[k + capasAntes][i + filasAntes][j + columnasAntes] = get(k, i, j);
        return crear(datos, mathContext, tolerancia);
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

    /** Combina elemento a elemento aplicando broadcasting sobre ejes unitarios. */
    public final T zipBroadcast(Matriz3DCompleja otra, BinaryOperator<BigComplex> funcion) {
        mismaPolitica(otra);
        Objects.requireNonNull(funcion, "La función no puede ser null");
        int[] a = {capas, filas, columnas};
        int[] b = {otra.capas(), otra.filas(), otra.columnas()};
        int[] salida = new int[3];
        for (int eje = 0; eje < 3; eje++) {
            if (a[eje] != b[eje] && a[eje] != 1 && b[eje] != 1)
                throw new DimensionMatrizException("Broadcast incompatible en el eje " + eje);
            salida[eje] = Math.max(a[eje], b[eje]);
        }
        BigComplex[][][] datos = new BigComplex[salida[0]][salida[1]][salida[2]];
        for (int k = 0; k < salida[0]; k++) for (int i = 0; i < salida[1]; i++) for (int j = 0; j < salida[2]; j++)
            datos[k][i][j] = redondear(funcion.apply(get(a[0] == 1 ? 0 : k, a[1] == 1 ? 0 : i, a[2] == 1 ? 0 : j),
                    otra.get(b[0] == 1 ? 0 : k, b[1] == 1 ? 0 : i, b[2] == 1 ? 0 : j)));
        return crear(datos, mathContext, tolerancia);
    }

    public final T suma(Matriz3DCompleja otra) { return zip(otra, (a, b) -> a.add(b, mathContext)); }
    public final T resta(Matriz3DCompleja otra) { return zip(otra, (a, b) -> a.subtract(b, mathContext)); }
    public final T sumaBroadcast(Matriz3DCompleja otra) { return zipBroadcast(otra, (a, b) -> a.add(b, mathContext)); }
    public final T restaBroadcast(Matriz3DCompleja otra) { return zipBroadcast(otra, (a, b) -> a.subtract(b, mathContext)); }
    public final T productoHadamardBroadcast(Matriz3DCompleja otra) {
        return zipBroadcast(otra, (a, b) -> a.multiply(b, mathContext));
    }
    public final T cocienteHadamardBroadcast(Matriz3DCompleja otra) {
        return zipBroadcast(otra, (a, b) -> {
            if (esCero(b)) throw new DivisionPorCeroException("Divisor broadcast nulo");
            return a.divide(b, mathContext);
        });
    }
    public final T opuesta() { return map(BigComplex::negate); }
    public final T conjugada() { return map(BigComplex::conjugate); }
    public final T parteReal() { return map(z -> BigComplex.valueOf(z.re)); }
    public final T parteImaginaria() { return map(z -> BigComplex.valueOf(z.im)); }
    public final T valorAbsoluto() { return map(z -> BigComplex.valueOf(z.abs(mathContext))); }
    public final T potenciaElementos(int exponente) {
        if (exponente == 0) return rellenar(BigComplex.ONE);
        if (exponente < 0) {
            for (BigComplex valor : flatten()) if (esCero(valor))
                throw new DivisionPorCeroException("Potencia negativa de un elemento nulo");
        }
        return map(z -> potenciaCompleja(z, exponente));
    }
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

    /**
     * t-product de Kilmer-Martin: producto matricial de tubos con convolución
     * circular sobre el eje de capas.
     */
    public final T productoT(Matriz3DCompleja otra) {
        mismaPolitica(otra);
        if (capas != otra.capas() || columnas != otra.filas())
            throw new DimensionMatrizException("t-product incompatible");
        BigComplex[][][] datos = nuevosDatos(capas, filas, otra.columnas(), BigComplex.ZERO);
        for (int k = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < otra.columnas(); j++) {
            BigComplex suma = BigComplex.ZERO;
            for (int p = 0; p < capas; p++) for (int n = 0; n < columnas; n++)
                suma = suma.add(get(p, i, n).multiply(otra.get(Math.floorMod(k - p, capas), n, j), mathContext), mathContext);
            datos[k][i][j] = suma.round(mathContext);
        }
        return crear(datos, mathContext, tolerancia);
    }

    /** Transformada discreta de Fourier unidimensional sobre las capas. */
    public final T transformadaFourierCapas() { return transformarFourierCapas(false); }
    public final T transformadaInversaFourierCapas() { return transformarFourierCapas(true); }

    private T transformarFourierCapas(boolean inversa) {
        BigComplex[][][] datos = nuevosDatos(capas, filas, columnas, BigComplex.ZERO);
        BigDecimal signo = inversa ? BigDecimal.ONE : BigDecimal.ONE.negate();
        BigDecimal dosPi = BigDecimalMath.pi(mathContext).multiply(BigDecimal.TWO, mathContext);
        BigDecimal escala = inversa ? BigDecimal.ONE.divide(BigDecimal.valueOf(capas), mathContext) : BigDecimal.ONE;
        for (int frecuencia = 0; frecuencia < capas; frecuencia++) for (int i = 0; i < filas; i++)
            for (int j = 0; j < columnas; j++) {
                BigComplex suma = BigComplex.ZERO;
                for (int k = 0; k < capas; k++) {
                    BigDecimal angulo = dosPi.multiply(BigDecimal.valueOf((long) frecuencia * k), mathContext)
                            .divide(BigDecimal.valueOf(capas), mathContext).multiply(signo, mathContext);
                    BigComplex raiz = BigComplex.valueOfPolar(BigDecimal.ONE, angulo, mathContext);
                    suma = suma.add(get(k, i, j).multiply(raiz, mathContext), mathContext);
                }
                datos[frecuencia][i][j] = suma.multiply(escala, mathContext).round(mathContext);
            }
        return crear(datos, mathContext, tolerancia);
    }

    /** Transpuesta conjugada en el álgebra t-product. */
    public final T adjuntaT() {
        BigComplex[][][] datos = new BigComplex[capas][columnas][filas];
        for (int k = 0; k < capas; k++) {
            int destino = k == 0 ? 0 : capas - k;
            for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
                datos[destino][j][i] = get(k, i, j).conjugate();
        }
        return crear(datos, mathContext, tolerancia);
    }

    /** Inversa respecto del t-product, calculada por capas en Fourier. */
    public final T inversaT() {
        if (filas != columnas) throw new DimensionMatrizException("La inversa T exige capas matriciales cuadradas");
        return transformadaFourierCapas().aplicarPorCapas(m -> m.inversa()).transformadaInversaFourierCapas();
    }

    /** t-SVD reducida obtenida mediante SVD matricial en el dominio de Fourier. */
    public final ResultadoTSVD<T> descomposicionTSVD() {
        T fourier = transformadaFourierCapas();
        List<P> us = new ArrayList<>(capas), sigmas = new ArrayList<>(capas), vs = new ArrayList<>(capas);
        for (int k = 0; k < capas; k++) {
            ResultadoSVD<P> svd = fourier.capa(k).descomposicionSVD();
            us.add(svd.u()); sigmas.add(svd.sigma()); vs.add(svd.v());
        }
        T u = ensamblarResultados(us).transformadaInversaFourierCapas();
        T sigma = ensamblarResultados(sigmas).transformadaInversaFourierCapas();
        T v = ensamblarResultados(vs).transformadaInversaFourierCapas();
        T reconstruida = u.productoT(sigma).productoT(v.adjuntaT());
        return new ResultadoTSVD<>(u, sigma, v, resta(reconstruida).normaFrobenius());
    }

    public final T reconstruirTSVD(ResultadoTSVD<? extends Matriz3DCompleja> resultado) {
        Objects.requireNonNull(resultado, "El resultado no puede ser null");
        T u = crear(copiarDatos(resultado.u()), mathContext, tolerancia);
        T sigma = crear(copiarDatos(resultado.sigma()), mathContext, tolerancia);
        T v = crear(copiarDatos(resultado.v()), mathContext, tolerancia);
        return u.productoT(sigma).productoT(v.adjuntaT());
    }

    public final boolean verificarTSVD(ResultadoTSVD<? extends Matriz3DCompleja> resultado, BigDecimal tol) {
        return aproximadamenteIgual(reconstruirTSVD(resultado), tol);
    }

    public final T potenciaT(int exponente) {
        if (filas != columnas) throw new DimensionMatrizException("La potencia T exige capas matriciales cuadradas");
        T base = exponente < 0 ? inversaT() : crear(copiaDatos(), mathContext, tolerancia);
        long n = Math.abs((long) exponente);
        T resultado = identidadTInterna(capas, filas);
        while (n > 0) {
            if ((n & 1L) != 0) resultado = resultado.productoT(base);
            n >>= 1;
            if (n > 0) base = base.productoT(base);
        }
        return resultado;
    }

    public final boolean esIdentidadT() {
        return filas == columnas && aproximadamenteIgual(identidadTInterna(capas, filas), tolerancia);
    }

    private T identidadTInterna(int numeroCapas, int orden) {
        BigComplex[][][] datos = nuevosDatos(numeroCapas, orden, orden, BigComplex.ZERO);
        for (int i = 0; i < orden; i++) datos[0][i][i] = BigComplex.ONE;
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

    public final T intercambiarEjes(int ejeA, int ejeB) {
        validarEje(ejeA); validarEje(ejeB);
        int[] p = {0, 1, 2};
        int temporal = p[ejeA]; p[ejeA] = p[ejeB]; p[ejeB] = temporal;
        return permutarEjes(p[0], p[1], p[2]);
    }

    /** Mueve un eje a una nueva posición, desplazando los demás. */
    public final T moverEje(int origen, int destino) {
        validarEje(origen); validarEje(destino);
        List<Integer> ejes = new ArrayList<>(List.of(0, 1, 2));
        int eje = ejes.remove(origen);
        ejes.add(destino, eje);
        return permutarEjes(ejes.get(0), ejes.get(1), ejes.get(2));
    }

    /** Invierte el orden de elementos en un eje. */
    public final T voltear(int eje) {
        validarEje(eje);
        int[] forma = {capas, filas, columnas};
        BigComplex[][][] datos = new BigComplex[capas][filas][columnas];
        for (int k = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++) {
            int[] p = {k, i, j}; p[eje] = forma[eje] - 1 - p[eje];
            datos[k][i][j] = get(p[0], p[1], p[2]);
        }
        return crear(datos, mathContext, tolerancia);
    }

    /** Desplazamiento circular sobre un eje; admite desplazamientos negativos. */
    public final T desplazarCiclicamente(int desplazamiento, int eje) {
        validarEje(eje);
        int[] forma = {capas, filas, columnas};
        int d = Math.floorMod(desplazamiento, forma[eje]);
        BigComplex[][][] datos = new BigComplex[capas][filas][columnas];
        for (int k = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++) {
            int[] salida = {k, i, j};
            int[] origen = salida.clone(); origen[eje] = Math.floorMod(salida[eje] - d, forma[eje]);
            datos[k][i][j] = get(origen[0], origen[1], origen[2]);
        }
        return crear(datos, mathContext, tolerancia);
    }

    /** Rota múltiplos de 90 grados en el plano definido por dos ejes. */
    public final T rotar90(int veces, int ejeA, int ejeB) {
        validarEje(ejeA); validarEje(ejeB);
        if (ejeA == ejeB) throw new IllegalArgumentException("Los ejes de rotación deben ser distintos");
        int giros = Math.floorMod(veces, 4);
        T resultado = crear(copiaDatos(), mathContext, tolerancia);
        for (int giro = 0; giro < giros; giro++) resultado = resultado.rotar90UnaVez(ejeA, ejeB);
        return resultado;
    }

    protected final T rotar90UnaVez(int ejeA, int ejeB) {
        int[] origenForma = {capas, filas, columnas};
        int[] salidaForma = origenForma.clone();
        salidaForma[ejeA] = origenForma[ejeB];
        salidaForma[ejeB] = origenForma[ejeA];
        BigComplex[][][] datos = new BigComplex[salidaForma[0]][salidaForma[1]][salidaForma[2]];
        for (int k = 0; k < salidaForma[0]; k++) for (int i = 0; i < salidaForma[1]; i++) for (int j = 0; j < salidaForma[2]; j++) {
            int[] salida = {k, i, j};
            int[] origen = salida.clone();
            origen[ejeA] = salida[ejeB];
            origen[ejeB] = origenForma[ejeB] - 1 - salida[ejeA];
            datos[k][i][j] = get(origen[0], origen[1], origen[2]);
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

    public final BigComplex productoElementos() {
        BigComplex producto = BigComplex.ONE;
        for (BigComplex valor : flatten()) producto = producto.multiply(valor, mathContext);
        return producto.round(mathContext);
    }

    /** Reduce un eje y devuelve una matriz sobre los otros dos ejes. */
    public final P reducirEje(int eje, BinaryOperator<BigComplex> reduccion) {
        validarEje(eje);
        Objects.requireNonNull(reduccion, "La reducción no puede ser null");
        int[][] formas = {{filas, columnas}, {capas, columnas}, {capas, filas}};
        BigComplex[][] datos = new BigComplex[formas[eje][0]][formas[eje][1]];
        for (int i = 0; i < datos.length; i++) for (int j = 0; j < datos[i].length; j++) {
            BigComplex[] fibra = fibra(eje, i, j);
            BigComplex acumulado = fibra[0];
            for (int n = 1; n < fibra.length; n++) acumulado = redondear(reduccion.apply(acumulado, fibra[n]));
            datos[i][j] = acumulado;
        }
        return crearPlano(datos, mathContext, tolerancia);
    }

    public final P sumaEje(int eje) { return reducirEje(eje, (a, b) -> a.add(b, mathContext)); }
    public final P productoEje(int eje) { return reducirEje(eje, (a, b) -> a.multiply(b, mathContext)); }
    public final P promedioEje(int eje) { return sumaEje(eje).dividirPorEscalar(BigDecimal.valueOf(tamanoEje(eje))); }

    public final P promedioPonderadoEje(int eje, BigDecimal[] pesos) {
        validarEje(eje);
        Objects.requireNonNull(pesos, "Los pesos no pueden ser null");
        if (pesos.length != tamanoEje(eje)) throw new DimensionMatrizException("Número de pesos incompatible");
        BigDecimal total = BigDecimal.ZERO;
        for (BigDecimal peso : pesos) total = total.add(Objects.requireNonNull(peso), mathContext);
        if (total.abs().compareTo(tolerancia) <= 0) throw new DivisionPorCeroException("La suma de pesos es cero");
        BigComplex[] complejos = new BigComplex[pesos.length];
        for (int n = 0; n < pesos.length; n++) complejos[n] = BigComplex.valueOf(pesos[n].divide(total, mathContext));
        return contraerModo(complejos, eje);
    }

    /** Producto interno hermítico de cada par de fibras sobre un eje. */
    public final P productoInternoEje(Matriz3DCompleja otra, int eje) {
        mismasDimensiones(otra); validarEje(eje);
        int[][] formas = {{filas, columnas}, {capas, columnas}, {capas, filas}};
        BigComplex[][] datos = new BigComplex[formas[eje][0]][formas[eje][1]];
        for (int i = 0; i < datos.length; i++) for (int j = 0; j < datos[i].length; j++) {
            BigComplex[] a = fibra(eje, i, j), b = fibraDe(otra, eje, i, j);
            BigComplex suma = BigComplex.ZERO;
            for (int n = 0; n < a.length; n++)
                suma = suma.add(a[n].conjugate().multiply(b[n], mathContext), mathContext);
            datos[i][j] = suma.round(mathContext);
        }
        return crearPlano(datos, mathContext, tolerancia);
    }

    public final P normaFrobeniusEje(int eje) {
        return productoInternoEje(this, eje).map(z -> BigComplex.valueOf(z.re.max(BigDecimal.ZERO).sqrt(mathContext)));
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

    /** Varianza poblacional compleja definida como la media de |x-media|². */
    public final BigDecimal varianza() {
        BigComplex media = promedio();
        BigDecimal suma = BigDecimal.ZERO;
        for (BigComplex valor : flatten()) {
            BigDecimal modulo = valor.subtract(media, mathContext).abs(mathContext);
            suma = suma.add(modulo.multiply(modulo, mathContext), mathContext);
        }
        return suma.divide(BigDecimal.valueOf(numeroElementos()), mathContext);
    }

    public final BigDecimal desviacionEstandar() { return varianza().sqrt(mathContext); }

    public final P varianzaEje(int eje) {
        validarEje(eje);
        int[][] formas = {{filas, columnas}, {capas, columnas}, {capas, filas}};
        BigComplex[][] datos = new BigComplex[formas[eje][0]][formas[eje][1]];
        for (int i = 0; i < datos.length; i++) for (int j = 0; j < datos[i].length; j++) {
            BigComplex[] valores = fibra(eje, i, j);
            BigComplex media = BigComplex.ZERO;
            for (BigComplex valor : valores) media = media.add(valor, mathContext);
            media = media.divide(BigDecimal.valueOf(valores.length), mathContext);
            BigDecimal suma = BigDecimal.ZERO;
            for (BigComplex valor : valores) {
                BigDecimal modulo = valor.subtract(media, mathContext).abs(mathContext);
                suma = suma.add(modulo.multiply(modulo, mathContext), mathContext);
            }
            datos[i][j] = BigComplex.valueOf(suma.divide(BigDecimal.valueOf(valores.length), mathContext));
        }
        return crearPlano(datos, mathContext, tolerancia);
    }

    public final int contarNoNulos() {
        int cantidad = 0;
        for (BigComplex valor : flatten()) if (!esCero(valor)) cantidad++;
        return cantidad;
    }

    public final BigComplex elementoMayorModulo() { int[] p = indiceMayorModulo(); return get(p[0], p[1], p[2]); }
    public final BigComplex elementoMenorModulo() { int[] p = indiceMenorModulo(); return get(p[0], p[1], p[2]); }
    public final int[] indiceMayorModulo() { return indiceExtremoModulo(true); }
    public final int[] indiceMenorModulo() { return indiceExtremoModulo(false); }

    private int[] indiceExtremoModulo(boolean mayor) {
        int[] mejor = {0, 0, 0};
        BigDecimal extremo = get(0, 0, 0).abs(mathContext);
        for (int k = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++) {
            BigDecimal actual = get(k, i, j).abs(mathContext);
            if (mayor ? actual.compareTo(extremo) > 0 : actual.compareTo(extremo) < 0) {
                extremo = actual; mejor = new int[]{k, i, j};
            }
        }
        return mejor;
    }

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

    /** Matricización modo-n; las columnas restantes usan orden lexicográfico. */
    public final P desplegar(int eje) {
        validarEje(eje);
        int[] forma = {capas, filas, columnas};
        BigComplex[][] datos = new BigComplex[forma[eje]][numeroElementos() / forma[eje]];
        for (int k = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++) {
            int fila = eje == 0 ? k : eje == 1 ? i : j;
            int columna = eje == 0 ? i * columnas + j : eje == 1 ? k * columnas + j : k * filas + i;
            datos[fila][columna] = get(k, i, j);
        }
        return crearPlano(datos, mathContext, tolerancia);
    }

    /** Producto tensor-matriz modo-n: el operando debe medir J × tamaño(eje). */
    public final T productoModo(MatrizCompleja matriz, int eje) {
        validarEje(eje);
        Objects.requireNonNull(matriz, "La matriz no puede ser null");
        if (!mathContext.equals(matriz.mathContext()) || tolerancia.compareTo(matriz.tolerancia()) != 0)
            throw new OperacionNoDefinidaException("No se mezclan políticas de precisión distintas");
        if (matriz.columnas() != tamanoEje(eje)) throw new DimensionMatrizException("Producto modo-n incompatible");
        int[] salida = {capas, filas, columnas}; salida[eje] = matriz.filas();
        BigComplex[][][] datos = nuevosDatos(salida[0], salida[1], salida[2], BigComplex.ZERO);
        for (int k = 0; k < salida[0]; k++) for (int i = 0; i < salida[1]; i++) for (int j = 0; j < salida[2]; j++) {
            int[] p = {k, i, j};
            BigComplex suma = BigComplex.ZERO;
            for (int n = 0; n < tamanoEje(eje); n++) {
                int salidaEje = p[eje]; p[eje] = n;
                suma = suma.add(matriz.get(salidaEje, n).multiply(get(p[0], p[1], p[2]), mathContext), mathContext);
                p[eje] = salidaEje;
            }
            datos[k][i][j] = suma.round(mathContext);
        }
        return crear(datos, mathContext, tolerancia);
    }

    public final T productoMultilineal(MatrizCompleja modo0, MatrizCompleja modo1, MatrizCompleja modo2) {
        return productoModo(modo0, 0).productoModo(modo1, 1).productoModo(modo2, 2);
    }

    /** Producto modo-n por vector; elimina el eje contraído. */
    public final P contraerModo(BigComplex[] vector, int eje) {
        validarEje(eje);
        Objects.requireNonNull(vector, "El vector no puede ser null");
        if (vector.length != tamanoEje(eje)) throw new DimensionMatrizException("Vector incompatible con el eje");
        int[][] formas = {{filas, columnas}, {capas, columnas}, {capas, filas}};
        BigComplex[][] datos = new BigComplex[formas[eje][0]][formas[eje][1]];
        for (int i = 0; i < datos.length; i++) for (int j = 0; j < datos[i].length; j++) {
            BigComplex[] valores = fibra(eje, i, j);
            BigComplex suma = BigComplex.ZERO;
            for (int n = 0; n < valores.length; n++)
                suma = suma.add(valores[n].multiply(Objects.requireNonNull(vector[n]), mathContext), mathContext);
            datos[i][j] = suma.round(mathContext);
        }
        return crearPlano(datos, mathContext, tolerancia);
    }

    /** Forma multilineal X ×0 a ×1 b ×2 c. */
    public final BigComplex formaMultilineal(BigComplex[] modo0, BigComplex[] modo1, BigComplex[] modo2) {
        validarLongitudVector(modo0, capas); validarLongitudVector(modo1, filas); validarLongitudVector(modo2, columnas);
        BigComplex suma = BigComplex.ZERO;
        for (int k = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            suma = suma.add(get(k, i, j).multiply(modo0[k], mathContext)
                    .multiply(modo1[i], mathContext).multiply(modo2[j], mathContext), mathContext);
        return suma.round(mathContext);
    }

    /** Rango multilineal: rangos de los tres desplegados. */
    public final int[] rangoMultilineal() {
        return new int[]{desplegar(0).rango(), desplegar(1).rango(), desplegar(2).rango()};
    }

    /** HOSVD completa mediante las SVD de los tres desplegados modo-n. */
    public final ResultadoHOSVD<T, P> descomposicionHOSVD() {
        int elementos = numeroElementos();
        return descomposicionHOSVDTruncada(Math.min(capas, elementos / capas),
                Math.min(filas, elementos / filas), Math.min(columnas, elementos / columnas));
    }

    /** Descomposición Tucker truncada a los rangos multilineales solicitados. */
    public final ResultadoHOSVD<T, P> descomposicionHOSVDTruncada(int rango0, int rango1, int rango2) {
        int[] rangos = {rango0, rango1, rango2};
        List<P> factores = factoresHOSVD(rangos);
        T nucleo = productoModo(factores.get(0).adjuntaHermitica(), 0)
                .productoModo(factores.get(1).adjuntaHermitica(), 1)
                .productoModo(factores.get(2).adjuntaHermitica(), 2);
        T reconstruida = nucleo.productoModo(factores.get(0), 0)
                .productoModo(factores.get(1), 1).productoModo(factores.get(2), 2);
        BigDecimal error = resta(reconstruida).normaFrobenius();
        return new ResultadoHOSVD<>(nucleo, factores.get(0), factores.get(1), factores.get(2), error);
    }

    public final T reconstruirHOSVD(ResultadoHOSVD<? extends Matriz3DCompleja, ? extends MatrizCompleja> resultado) {
        Objects.requireNonNull(resultado, "El resultado no puede ser null");
        Matriz3DCompleja nucleo = resultado.nucleo();
        T base = crear(copiarDatos(nucleo), mathContext, tolerancia);
        return base.productoModo(resultado.factorCapas(), 0)
                .productoModo(resultado.factorFilas(), 1).productoModo(resultado.factorColumnas(), 2);
    }

    public final boolean verificarHOSVD(ResultadoHOSVD<? extends Matriz3DCompleja, ? extends MatrizCompleja> resultado,
                                         BigDecimal tol) {
        return aproximadamenteIgual(reconstruirHOSVD(resultado), tol);
    }

    private List<P> factoresHOSVD(int[] rangos) {
        List<P> factores = new ArrayList<>(3);
        for (int eje = 0; eje < 3; eje++) {
            P desplegada = desplegar(eje);
            int maximo = Math.min(desplegada.filas(), desplegada.columnas());
            if (rangos[eje] <= 0 || rangos[eje] > maximo)
                throw new DimensionMatrizException("Rango HOSVD inválido para el eje " + eje);
            try {
                P gramIzquierda = desplegada.producto(desplegada.adjuntaHermitica());
                ResultadoDiagonalizacion<P> espectro = gramIzquierda.diagonalizacionUnitaria();
                List<Integer> orden = new ArrayList<>(gramIzquierda.filas());
                for (int n = 0; n < gramIzquierda.filas(); n++) orden.add(n);
                orden.sort((a, b) -> espectro.d().get(b, b).re.compareTo(espectro.d().get(a, a).re));
                BigComplex[][] factor = new BigComplex[gramIzquierda.filas()][rangos[eje]];
                for (int columna = 0; columna < rangos[eje]; columna++)
                    for (int fila = 0; fila < gramIzquierda.filas(); fila++)
                        factor[fila][columna] = espectro.p().get(fila, orden.get(columna));
                factores.add(crearPlano(factor, mathContext, tolerancia));
            } catch (MatrizException excepcionNumerica) {
                factores.add(baseCanonicaPorEnergia(desplegada, rangos[eje]));
            }
        }
        return List.copyOf(factores);
    }

    /** Respaldo ortonormal determinista para desplegados numéricamente degenerados. */
    private P baseCanonicaPorEnergia(P desplegada, int rango) {
        List<Integer> indices = new ArrayList<>(desplegada.filas());
        for (int i = 0; i < desplegada.filas(); i++) indices.add(i);
        indices.sort((a, b) -> energiaFila(desplegada, b).compareTo(energiaFila(desplegada, a)));
        BigComplex[][] base = new BigComplex[desplegada.filas()][rango];
        for (BigComplex[] fila : base) Arrays.fill(fila, BigComplex.ZERO);
        for (int columna = 0; columna < rango; columna++) base[indices.get(columna)][columna] = BigComplex.ONE;
        return crearPlano(base, mathContext, tolerancia);
    }

    private BigDecimal energiaFila(P matriz, int fila) {
        BigDecimal suma = BigDecimal.ZERO;
        for (int j = 0; j < matriz.columnas(); j++) {
            BigDecimal modulo = matriz.get(fila, j).abs(mathContext);
            suma = suma.add(modulo.multiply(modulo, mathContext), mathContext);
        }
        return suma;
    }

    /** Traza entre dos ejes de igual longitud; conserva el eje restante. */
    public final BigComplex[] traza(int ejeA, int ejeB) {
        validarEje(ejeA); validarEje(ejeB);
        if (ejeA == ejeB) throw new IllegalArgumentException("Los ejes de traza deben ser distintos");
        if (tamanoEje(ejeA) != tamanoEje(ejeB)) throw new DimensionMatrizException("Los ejes de traza deben medir igual");
        int restante = 3 - ejeA - ejeB;
        BigComplex[] salida = new BigComplex[tamanoEje(restante)];
        for (int r = 0; r < salida.length; r++) {
            BigComplex suma = BigComplex.ZERO;
            for (int n = 0; n < tamanoEje(ejeA); n++) {
                int[] p = new int[3]; p[ejeA] = n; p[ejeB] = n; p[restante] = r;
                suma = suma.add(get(p[0], p[1], p[2]), mathContext);
            }
            salida[r] = suma.round(mathContext);
        }
        return salida;
    }

    public final BigComplex[] diagonalPrincipal() {
        if (capas != filas || filas != columnas) throw new DimensionMatrizException("La diagonal principal exige tensor cúbico");
        BigComplex[] diagonal = new BigComplex[capas];
        for (int n = 0; n < capas; n++) diagonal[n] = get(n, n, n);
        return diagonal;
    }

    public final T productoKronecker(Matriz3DCompleja otra) {
        mismaPolitica(otra);
        int nc = Math.multiplyExact(capas, otra.capas());
        int nf = Math.multiplyExact(filas, otra.filas());
        int nn = Math.multiplyExact(columnas, otra.columnas());
        BigComplex[][][] datos = new BigComplex[nc][nf][nn];
        for (int k = 0; k < nc; k++) for (int i = 0; i < nf; i++) for (int j = 0; j < nn; j++)
            datos[k][i][j] = get(k / otra.capas(), i / otra.filas(), j / otra.columnas())
                    .multiply(otra.get(k % otra.capas(), i % otra.filas(), j % otra.columnas()), mathContext).round(mathContext);
        return crear(datos, mathContext, tolerancia);
    }

    /** Producto vectorial de fibras; el eje elegido debe medir tres. */
    public final T productoVectorial(Matriz3DCompleja otra, int eje) {
        mismasDimensiones(otra); validarEje(eje);
        if (tamanoEje(eje) != 3) throw new DimensionMatrizException("El producto vectorial exige longitud tres");
        BigComplex[][][] datos = nuevosDatos(capas, filas, columnas, BigComplex.ZERO);
        int[][] formas = {{filas, columnas}, {capas, columnas}, {capas, filas}};
        for (int a = 0; a < formas[eje][0]; a++) for (int b = 0; b < formas[eje][1]; b++) {
            BigComplex[] x = fibra(eje, a, b);
            BigComplex[] y = fibraDe(otra, eje, a, b);
            BigComplex[] z = {
                    x[1].multiply(y[2], mathContext).subtract(x[2].multiply(y[1], mathContext), mathContext),
                    x[2].multiply(y[0], mathContext).subtract(x[0].multiply(y[2], mathContext), mathContext),
                    x[0].multiply(y[1], mathContext).subtract(x[1].multiply(y[0], mathContext), mathContext)};
            for (int n = 0; n < 3; n++) { int[] p = coordenadasFibra(eje, n, a, b); datos[p[0]][p[1]][p[2]] = z[n].round(mathContext); }
        }
        return crear(datos, mathContext, tolerancia);
    }

    /** Doble contracción tensorial; devuelve los dos ejes no contraídos. */
    public final P contraerDoble(Matriz3DCompleja otra, int ejeA1, int ejeA2, int ejeB1, int ejeB2) {
        mismaPolitica(otra);
        validarDosEjes(ejeA1, ejeA2); validarDosEjes(ejeB1, ejeB2);
        int[] a = {capas, filas, columnas}; int[] b = {otra.capas(), otra.filas(), otra.columnas()};
        if (a[ejeA1] != b[ejeB1] || a[ejeA2] != b[ejeB2])
            throw new DimensionMatrizException("Ejes incompatibles para la doble contracción");
        int libreA = 3 - ejeA1 - ejeA2; int libreB = 3 - ejeB1 - ejeB2;
        BigComplex[][] datos = new BigComplex[a[libreA]][b[libreB]];
        for (int i = 0; i < datos.length; i++) for (int j = 0; j < datos[i].length; j++) {
            BigComplex suma = BigComplex.ZERO;
            for (int p = 0; p < a[ejeA1]; p++) for (int q = 0; q < a[ejeA2]; q++) {
                int[] ca = new int[3]; ca[libreA] = i; ca[ejeA1] = p; ca[ejeA2] = q;
                int[] cb = new int[3]; cb[libreB] = j; cb[ejeB1] = p; cb[ejeB2] = q;
                suma = suma.add(get(ca[0], ca[1], ca[2]).multiply(otra.get(cb[0], cb[1], cb[2]), mathContext), mathContext);
            }
            datos[i][j] = suma.round(mathContext);
        }
        return crearPlano(datos, mathContext, tolerancia);
    }

    /** Producto exterior sin contracciones; la salida tiene orden seis. */
    public final BigComplex[][][][][][] productoExterior(Matriz3DCompleja otra) {
        mismaPolitica(otra);
        BigComplex[][][][][][] datos = new BigComplex[capas][filas][columnas]
                [otra.capas()][otra.filas()][otra.columnas()];
        for (int k = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            for (int p = 0; p < otra.capas(); p++) for (int q = 0; q < otra.filas(); q++) for (int r = 0; r < otra.columnas(); r++)
                datos[k][i][j][p][q][r] = get(k, i, j).multiply(otra.get(p, q, r), mathContext).round(mathContext);
        return datos;
    }

    /** Convolución lineal discreta tridimensional directa. */
    public final T convolucion(Matriz3DCompleja nucleo, ModoConvolucion modo) {
        mismaPolitica(nucleo);
        Objects.requireNonNull(modo, "El modo no puede ser null");
        if (modo == ModoConvolucion.VALIDA && (capas < nucleo.capas() || filas < nucleo.filas()
                || columnas < nucleo.columnas()))
            throw new DimensionMatrizException("La convolución válida exige que el tensor contenga al núcleo");
        int fc = Math.addExact(capas, nucleo.capas()) - 1;
        int ff = Math.addExact(filas, nucleo.filas()) - 1;
        int fn = Math.addExact(columnas, nucleo.columnas()) - 1;
        BigComplex[][][] completa = nuevosDatos(fc, ff, fn, BigComplex.ZERO);
        for (int ok = 0; ok < fc; ok++) for (int oi = 0; oi < ff; oi++) for (int oj = 0; oj < fn; oj++) {
            BigComplex suma = BigComplex.ZERO;
            for (int kk = 0; kk < nucleo.capas(); kk++) for (int ki = 0; ki < nucleo.filas(); ki++)
                for (int kj = 0; kj < nucleo.columnas(); kj++) {
                    int ak = ok - kk, ai = oi - ki, aj = oj - kj;
                    if (ak >= 0 && ak < capas && ai >= 0 && ai < filas && aj >= 0 && aj < columnas)
                        suma = suma.add(get(ak, ai, aj).multiply(nucleo.get(kk, ki, kj), mathContext), mathContext);
                }
            completa[ok][oi][oj] = suma.round(mathContext);
        }
        T total = crear(completa, mathContext, tolerancia);
        if (modo == ModoConvolucion.COMPLETA) return total;
        if (modo == ModoConvolucion.MISMA) {
            int ck = (nucleo.capas() - 1) / 2;
            int ci = (nucleo.filas() - 1) / 2;
            int cj = (nucleo.columnas() - 1) / 2;
            return total.submatriz(ck, ck + capas, ci, ci + filas, cj, cj + columnas);
        }
        return total.submatriz(nucleo.capas() - 1, capas,
                nucleo.filas() - 1, filas, nucleo.columnas() - 1, columnas);
    }

    /** Correlación cruzada 3D, conjugando el segundo operando complejo. */
    public final T correlacion(Matriz3DCompleja nucleo, ModoConvolucion modo) {
        mismaPolitica(nucleo);
        BigComplex[][][] reflejado = new BigComplex[nucleo.capas()][nucleo.filas()][nucleo.columnas()];
        for (int k = 0; k < nucleo.capas(); k++) for (int i = 0; i < nucleo.filas(); i++)
            for (int j = 0; j < nucleo.columnas(); j++)
                reflejado[k][i][j] = nucleo.get(nucleo.capas() - 1 - k, nucleo.filas() - 1 - i,
                        nucleo.columnas() - 1 - j).conjugate();
        return convolucion(crear(reflejado, mathContext, tolerancia), modo);
    }

    public final boolean esCubica() { return capas == filas && filas == columnas; }

    public final boolean esDiagonal() {
        for (int k = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            if (!(k == i && i == j) && !esCero(get(k, i, j))) return false;
        return true;
    }

    /** Simetría bajo cualquier permutación de ejes; exige tensor cúbico. */
    public final boolean esSimetrica() {
        if (!esCubica()) return false;
        return igualExacta(permutarEjes(0, 2, 1)) && igualExacta(permutarEjes(1, 0, 2));
    }

    /** Antisimetría bajo el intercambio de cualquier par de ejes. */
    public final boolean esAntisimetrica() {
        if (!esCubica()) return false;
        return opuesta().igualExacta(permutarEjes(0, 2, 1))
                && opuesta().igualExacta(permutarEjes(1, 0, 2));
    }

    /** Promedio sobre las seis permutaciones de ejes. */
    public final T simetrizar() {
        if (!esCubica()) throw new DimensionMatrizException("La simetrización exige tensor cúbico");
        int[][] permutaciones = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
        T suma = rellenar(BigComplex.ZERO);
        for (int[] p : permutaciones) suma = suma.suma(permutarEjes(p[0], p[1], p[2]));
        return suma.dividirPorEscalar(BigDecimal.valueOf(6));
    }

    /** Suma alternada sobre las seis permutaciones de ejes. */
    public final T antisimetrizar() {
        if (!esCubica()) throw new DimensionMatrizException("La antisimetrización exige tensor cúbico");
        int[][] p = {{0,1,2},{1,2,0},{2,0,1},{0,2,1},{2,1,0},{1,0,2}};
        T suma = rellenar(BigComplex.ZERO);
        for (int n = 0; n < p.length; n++) {
            T termino = permutarEjes(p[n][0], p[n][1], p[n][2]);
            suma = n < 3 ? suma.suma(termino) : suma.resta(termino);
        }
        return suma.dividirPorEscalar(BigDecimal.valueOf(6));
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

    protected static BigComplex[][][] datosDesdeVector(BigComplex[] valores, int capas, int filas, int columnas) {
        Objects.requireNonNull(valores, "Los valores no pueden ser null");
        int total;
        try { total = Math.multiplyExact(capas, Math.multiplyExact(filas, columnas)); }
        catch (ArithmeticException e) { throw new DimensionMatrizException("Dimensiones demasiado grandes"); }
        if (capas <= 0 || filas <= 0 || columnas <= 0 || valores.length != total)
            throw new DimensionMatrizException("El vector no coincide con las dimensiones 3D");
        BigComplex[][][] datos = new BigComplex[capas][filas][columnas];
        for (int k = 0, n = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            datos[k][i][j] = Objects.requireNonNull(valores[n++], "Un elemento no puede ser null");
        return datos;
    }

    protected static BigComplex[][][] parsearSerializado(String texto) {
        Objects.requireNonNull(texto, "El texto no puede ser null");
        int separador = texto.indexOf(':');
        if (separador <= 0) throw new IllegalArgumentException("Formato 3D inválido");
        String[] forma = texto.substring(0, separador).split("x", -1);
        if (forma.length != 3) throw new IllegalArgumentException("Formato de dimensiones inválido");
        int capas, filas, columnas;
        try {
            capas = Integer.parseInt(forma[0]); filas = Integer.parseInt(forma[1]); columnas = Integer.parseInt(forma[2]);
        } catch (NumberFormatException e) { throw new IllegalArgumentException("Dimensiones inválidas", e); }
        String cuerpo = texto.substring(separador + 1);
        String[] elementos = cuerpo.isEmpty() ? new String[0] : cuerpo.split(";", -1);
        BigComplex[] valores = new BigComplex[elementos.length];
        try {
            for (int n = 0; n < elementos.length; n++) {
                String[] partes = elementos[n].split(",", -1);
                if (partes.length != 2) throw new IllegalArgumentException("Elemento complejo inválido");
                valores[n] = BigComplex.valueOf(new BigDecimal(partes[0]), new BigDecimal(partes[1]));
            }
        } catch (NumberFormatException e) { throw new IllegalArgumentException("Número decimal inválido", e); }
        return datosDesdeVector(valores, capas, filas, columnas);
    }

    private P capaDe(Matriz3DCompleja origen, int capa) {
        BigComplex[][] datos = new BigComplex[origen.filas()][origen.columnas()];
        for (int i = 0; i < origen.filas(); i++)
            for (int j = 0; j < origen.columnas(); j++) datos[i][j] = origen.get(capa, i, j);
        return crearPlano(datos, mathContext, tolerancia);
    }

    private static BigComplex[][][] copiarDatos(Matriz3DCompleja origen) {
        BigComplex[][][] datos = new BigComplex[origen.capas()][origen.filas()][origen.columnas()];
        for (int k = 0; k < origen.capas(); k++) for (int i = 0; i < origen.filas(); i++)
            for (int j = 0; j < origen.columnas(); j++) datos[k][i][j] = origen.get(k, i, j);
        return datos;
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

    private T intervaloEje(int eje, int desde, int hasta) {
        return switch (eje) {
            case 0 -> submatriz(desde, hasta, 0, filas, 0, columnas);
            case 1 -> submatriz(0, capas, desde, hasta, 0, columnas);
            default -> submatriz(0, capas, 0, filas, desde, hasta);
        };
    }

    private static int[] coordenadasFibra(int eje, int variable, int indiceA, int indiceB) {
        return switch (eje) {
            case 0 -> new int[]{variable, indiceA, indiceB};
            case 1 -> new int[]{indiceA, variable, indiceB};
            case 2 -> new int[]{indiceA, indiceB, variable};
            default -> throw new IllegalArgumentException("El eje debe ser 0, 1 o 2");
        };
    }

    private static int[] coordenadasPlano(int eje, int capa, int fila, int columna) {
        return switch (eje) {
            case 0 -> new int[]{fila, columna};
            case 1 -> new int[]{capa, columna};
            case 2 -> new int[]{capa, fila};
            default -> throw new IllegalArgumentException("El eje debe ser 0, 1 o 2");
        };
    }

    private BigComplex[] fibraDe(Matriz3DCompleja origen, int eje, int indiceA, int indiceB) {
        BigComplex[] salida = new BigComplex[tamanoEje(eje)];
        for (int n = 0; n < salida.length; n++) {
            int[] p = coordenadasFibra(eje, n, indiceA, indiceB);
            salida[n] = origen.get(p[0], p[1], p[2]);
        }
        return salida;
    }

    private BigComplex potenciaCompleja(BigComplex base, int exponente) {
        long n = Math.abs((long) exponente);
        BigComplex resultado = BigComplex.ONE;
        BigComplex factor = base;
        while (n > 0) {
            if ((n & 1L) != 0) resultado = resultado.multiply(factor, mathContext);
            n >>= 1;
            if (n > 0) factor = factor.multiply(factor, mathContext);
        }
        return exponente < 0 ? BigComplex.ONE.divide(resultado, mathContext) : resultado;
    }

    private static void validarLongitudVector(BigComplex[] vector, int longitud) {
        Objects.requireNonNull(vector, "El vector no puede ser null");
        if (vector.length != longitud) throw new DimensionMatrizException("Longitud de vector incompatible");
        for (BigComplex valor : vector) Objects.requireNonNull(valor, "Un elemento del vector no puede ser null");
    }

    private static void validarDosEjes(int ejeA, int ejeB) {
        validarEje(ejeA); validarEje(ejeB);
        if (ejeA == ejeB) throw new IllegalArgumentException("Los ejes contraídos deben ser distintos");
    }

    private static void validarEje(int eje) {
        if (eje < 0 || eje > 2) throw new IllegalArgumentException("El eje debe ser 0, 1 o 2");
    }

    private static void validarPermutacion(int[] p) {
        if (Arrays.stream(p).anyMatch(x -> x < 0 || x > 2) || Arrays.stream(p).distinct().count() != 3)
            throw new IllegalArgumentException("Los ejes deben ser una permutación de 0, 1 y 2");
    }
}
