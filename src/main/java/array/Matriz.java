package array;

import ch.obermuhlner.math.big.BigComplex;
import comun.AbstractMatriz;
import comun.Complejos;
import comun.MatrizCompleja;
import comun.excepciones.DimensionMatrizException;
import comun.excepciones.OperacionNoDefinidaException;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Arrays;
import java.util.Objects;
import java.util.Random;

/**
 * Matriz inmutable almacenada realmente en {@code BigComplex[][]}. Los
 * constructores y getters copian profundamente los contenedores, BigComplex es
 * inmutable, y los algoritmos heredados devuelven siempre una nueva Matriz.
 * Se usa DECIMAL128 y tolerancia 1E-28 salvo configuración explícita.
 */
public final class Matriz extends AbstractMatriz<Matriz> {
    private final BigComplex[][] matriz;

    public Matriz(int filas, int columnas, boolean diagonal) {
        this(datosDimensiones(filas, columnas, diagonal ? BigComplex.ONE : BigComplex.ZERO, diagonal),
                Complejos.CONTEXTO_DEFECTO, Complejos.TOLERANCIA_DEFECTO);
    }
    public Matriz(double[][] array) { this(convertir(array), Complejos.CONTEXTO_DEFECTO, Complejos.TOLERANCIA_DEFECTO); }
    public Matriz(Double[][] array) { this(convertir(array), Complejos.CONTEXTO_DEFECTO, Complejos.TOLERANCIA_DEFECTO); }
    public Matriz(BigDecimal[][] array) { this(convertir(array), Complejos.CONTEXTO_DEFECTO, Complejos.TOLERANCIA_DEFECTO); }
    public Matriz(BigComplex[][] array) { this(array, Complejos.CONTEXTO_DEFECTO, Complejos.TOLERANCIA_DEFECTO); }
    public Matriz(String[][] array) { this(convertir(array), Complejos.CONTEXTO_DEFECTO, Complejos.TOLERANCIA_DEFECTO); }
    public Matriz(Matriz m) { this(Objects.requireNonNull(m).matriz, m.mathContext(), m.tolerancia()); }
    public Matriz(int filas, int columnas) { this(filas, columnas, false); }
    public Matriz(int filas, int columnas, BigComplex valorInicial) {
        this(datosDimensiones(filas, columnas, valorInicial, false), Complejos.CONTEXTO_DEFECTO, Complejos.TOLERANCIA_DEFECTO);
    }
    public Matriz(int filas, int columnas, MathContext mathContext) {
        this(datosDimensiones(filas, columnas, BigComplex.ZERO, false), mathContext,
                BigDecimal.ONE.scaleByPowerOfTen(-Math.max(1, mathContext.getPrecision() - 6)));
    }
    public Matriz(MatrizCompleja m) {
        this(m, m.mathContext(), m.tolerancia());
    }
    public Matriz(MatrizCompleja m, MathContext contexto, BigDecimal tolerancia) {
        this(copiar(m), contexto, tolerancia);
    }
    public Matriz(BigComplex[][] array, MathContext contexto, BigDecimal tolerancia) {
        this(validarYCopiar(array, contexto), contexto, tolerancia, true);
    }
    private Matriz(BigComplex[][] datos, MathContext contexto, BigDecimal tolerancia, boolean validados) {
        super(datos.length, datos[0].length, contexto, tolerancia);
        matriz = validados ? copiarEstructura(datos) : validarYCopiar(datos, contexto);
    }

    @Override public BigComplex get(int fila, int columna) {
        validarIndice(fila, columna);
        return matriz[fila][columna];
    }
    @Override protected Matriz crear(BigComplex[][] datos, MathContext contexto, BigDecimal tol) {
        return new Matriz(datos, contexto, tol);
    }

    public static Matriz cero(int filas, int columnas) { return new Matriz(filas, columnas); }
    public static Matriz identidad(int n) { return new Matriz(n, n, true); }
    public static Matriz identidad(int filas, int columnas) { return new Matriz(filas, columnas, true); }
    public static Matriz diagonal(BigComplex... elementos) {
        if (elementos == null || elementos.length == 0) throw new DimensionMatrizException("Diagonal vacía");
        BigComplex[][] d = datosDimensiones(elementos.length, elementos.length, BigComplex.ZERO, false);
        for (int i = 0; i < elementos.length; i++) d[i][i] = Objects.requireNonNull(elementos[i]);
        return new Matriz(d);
    }
    public static Matriz escalar(int n, BigComplex valor) {
        BigComplex[] d = new BigComplex[n]; Arrays.fill(d, Objects.requireNonNull(valor)); return diagonal(d);
    }
    public static Matriz fila(BigComplex... valores) { return new Matriz(new BigComplex[][]{valores}); }
    public static Matriz columna(BigComplex... valores) {
        if (valores == null || valores.length == 0) throw new DimensionMatrizException("Columna vacía");
        BigComplex[][] d = new BigComplex[valores.length][1];
        for (int i = 0; i < valores.length; i++) d[i][0] = valores[i];
        return new Matriz(d);
    }
    public static Matriz desdeFilas(BigComplex[]... filas) { return new Matriz(filas); }
    public static Matriz desdeColumnas(BigComplex[]... columnas) {
        if (columnas == null || columnas.length == 0 || columnas[0] == null)
            throw new DimensionMatrizException("Columnas vacías");
        BigComplex[][] d = new BigComplex[columnas[0].length][columnas.length];
        for (int j = 0; j < columnas.length; j++) {
            if (columnas[j] == null || columnas[j].length != d.length) throw new DimensionMatrizException("Columnas irregulares");
            for (int i = 0; i < d.length; i++) d[i][j] = columnas[j][i];
        }
        return new Matriz(d);
    }
    public static Matriz aleatoria(int filas, int columnas, long semilla) {
        Random r = new Random(semilla); BigComplex[][] d = new BigComplex[filas][columnas];
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            d[i][j] = BigComplex.valueOf(BigDecimal.valueOf(r.nextDouble() * 2 - 1));
        return new Matriz(d);
    }
    public static Matriz aleatoriaHermitica(int n, long semilla) {
        Matriz a = aleatoria(n, n, semilla); return a.suma(a.adjuntaHermitica()).dividirPorEscalar(BigDecimal.TWO);
    }
    public static Matriz aleatoriaUnitaria(int n, long semilla) {
        return aleatoria(n, n, semilla).descomposicionQRGramSchmidtModificado().q();
    }
    public static Matriz aleatoriaDefinidaPositiva(int n, long semilla) {
        Matriz a = aleatoria(n, n, semilla); return a.adjuntaHermitica().producto(a).suma(identidad(n));
    }
    public static Matriz aleatoriaDispersa(int filas, int columnas, BigDecimal densidad, long semilla) {
        if (densidad.signum() < 0 || densidad.compareTo(BigDecimal.ONE) > 0) throw new IllegalArgumentException("Densidad fuera de [0,1]");
        Random r = new Random(semilla); BigComplex[][] d = datosDimensiones(filas, columnas, BigComplex.ZERO, false);
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            if (r.nextDouble() < densidad.doubleValue()) d[i][j] = Complejos.de(r.nextDouble());
        return new Matriz(d);
    }
    public static Matriz matrizDePermutacion(int... permutacion) {
        if (permutacion == null || Arrays.stream(permutacion).distinct().count() != permutacion.length
                || Arrays.stream(permutacion).anyMatch(x -> x < 0 || x >= permutacion.length))
            throw new DimensionMatrizException("Permutación inválida");
        BigComplex[][] d = datosDimensiones(permutacion.length, permutacion.length, BigComplex.ZERO, false);
        for (int i = 0; i < permutacion.length; i++) d[i][permutacion[i]] = BigComplex.ONE;
        return new Matriz(d);
    }
    public static Matriz matrizDeRotacion2D(BigDecimal angulo) {
        MathContext mc = Complejos.CONTEXTO_DEFECTO;
        BigDecimal c = new BigDecimal(Math.cos(angulo.doubleValue()), mc);
        BigDecimal s = new BigDecimal(Math.sin(angulo.doubleValue()), mc);
        return new Matriz(new BigDecimal[][]{{c, s.negate()}, {s, c}});
    }
    public static Matriz matrizDeRotacion3D(BigDecimal angulo, char eje) {
        Matriz r = matrizDeRotacion2D(angulo);
        BigComplex c = r.get(0, 0), s = r.get(1, 0);
        return switch (Character.toLowerCase(eje)) {
            case 'x' -> new Matriz(new BigComplex[][]{{BigComplex.ONE, BigComplex.ZERO, BigComplex.ZERO},
                    {BigComplex.ZERO, c, s.negate()}, {BigComplex.ZERO, s, c}});
            case 'y' -> new Matriz(new BigComplex[][]{{c, BigComplex.ZERO, s}, {BigComplex.ZERO, BigComplex.ONE, BigComplex.ZERO},
                    {s.negate(), BigComplex.ZERO, c}});
            case 'z' -> new Matriz(new BigComplex[][]{{c, s.negate(), BigComplex.ZERO}, {s, c, BigComplex.ZERO},
                    {BigComplex.ZERO, BigComplex.ZERO, BigComplex.ONE}});
            default -> throw new IllegalArgumentException("Eje debe ser x, y o z");
        };
    }
    public static Matriz matrizDeReflexion(BigComplex[] normal) {
        Matriz n = columna(normal).normalizarColumnas();
        return identidad(normal.length).resta(n.producto(n.adjuntaHermitica()).productoEscalar(BigDecimal.TWO));
    }
    public static Matriz matrizDeProyeccion(BigComplex[] direccion) {
        Matriz n = columna(direccion).normalizarColumnas(); return n.producto(n.adjuntaHermitica());
    }
    public static Matriz householder(BigComplex[] vector) { return matrizDeReflexion(vector); }
    public static Matriz givens(int n, int i, int j, BigDecimal angulo) {
        Matriz g = identidad(n), r = matrizDeRotacion2D(angulo);
        return g.conElemento(i, i, r.get(0, 0)).conElemento(i, j, r.get(0, 1))
                .conElemento(j, i, r.get(1, 0)).conElemento(j, j, r.get(1, 1));
    }
    public static Matriz toeplitz(BigComplex[] primeraColumna, BigComplex[] primeraFila) {
        if (primeraColumna == null || primeraFila == null || primeraColumna.length == 0 || primeraFila.length == 0
                || !Objects.equals(primeraColumna[0], primeraFila[0])) throw new DimensionMatrizException("Datos Toeplitz inválidos");
        BigComplex[][] d = new BigComplex[primeraColumna.length][primeraFila.length];
        for (int i = 0; i < d.length; i++) for (int j = 0; j < d[i].length; j++)
            d[i][j] = j >= i ? primeraFila[j - i] : primeraColumna[i - j];
        return new Matriz(d);
    }
    public static Matriz hankel(BigComplex[] primeraFila, BigComplex[] ultimaColumna) {
        if (primeraFila == null || ultimaColumna == null || primeraFila.length == 0 || ultimaColumna.length == 0
                || !Objects.equals(primeraFila[primeraFila.length - 1], ultimaColumna[0]))
            throw new DimensionMatrizException("Datos Hankel inválidos");
        BigComplex[] sec = new BigComplex[primeraFila.length + ultimaColumna.length - 1];
        System.arraycopy(primeraFila, 0, sec, 0, primeraFila.length);
        System.arraycopy(ultimaColumna, 1, sec, primeraFila.length, ultimaColumna.length - 1);
        BigComplex[][] d = new BigComplex[ultimaColumna.length][primeraFila.length];
        for (int i = 0; i < d.length; i++) for (int j = 0; j < d[i].length; j++) d[i][j] = sec[i + j];
        return new Matriz(d);
    }
    public static Matriz circulante(BigComplex... primeraFila) {
        int n = primeraFila.length; BigComplex[][] d = new BigComplex[n][n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) d[i][j] = primeraFila[(j - i + n) % n];
        return new Matriz(d);
    }
    public static Matriz vandermonde(BigComplex... nodos) {
        BigComplex[][] d = new BigComplex[nodos.length][nodos.length];
        for (int i = 0; i < nodos.length; i++) {
            d[i][0] = BigComplex.ONE;
            for (int j = 1; j < nodos.length; j++) d[i][j] = d[i][j - 1].multiply(nodos[i], Complejos.CONTEXTO_DEFECTO);
        }
        return new Matriz(d);
    }
    public static Matriz hilbert(int n) {
        BigDecimal[][] d = new BigDecimal[n][n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
            d[i][j] = BigDecimal.ONE.divide(BigDecimal.valueOf(i + j + 1L), Complejos.CONTEXTO_DEFECTO);
        return new Matriz(d);
    }
    public static Matriz fourierDiscreta(int n) {
        BigComplex[][] d = new BigComplex[n][n]; MathContext mc = Complejos.CONTEXTO_DEFECTO;
        BigDecimal escala = BigDecimal.ONE.divide(BigDecimal.valueOf(n).sqrt(mc), mc);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
            d[i][j] = BigComplex.valueOfPolar(escala, new BigDecimal(-2 * Math.PI * i * j / n, mc), mc);
        return new Matriz(d);
    }
    public static Matriz hadamard(int n) {
        if (n <= 0 || (n & (n - 1)) != 0) throw new OperacionNoDefinidaException("Sylvester exige potencia de dos");
        BigComplex[][] d = {{BigComplex.ONE}};
        while (d.length < n) {
            int m = d.length; BigComplex[][] x = new BigComplex[2 * m][2 * m];
            for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) {
                x[i][j] = d[i][j]; x[i][j + m] = d[i][j]; x[i + m][j] = d[i][j]; x[i + m][j + m] = d[i][j].negate();
            }
            d = x;
        }
        return new Matriz(d);
    }
    public static Matriz triangularSuperior(BigComplex[][] datos) {
        Matriz m = new Matriz(datos);
        if (!m.esTriangularSuperior()) throw new DimensionMatrizException("Los datos no son triangulares superiores");
        return m;
    }
    public static Matriz triangularInferior(BigComplex[][] datos) {
        Matriz m = new Matriz(datos);
        if (!m.esTriangularInferior()) throw new DimensionMatrizException("Los datos no son triangulares inferiores");
        return m;
    }
    public static Matriz bidiagonalSuperior(BigComplex[] diagonal, BigComplex[] superior) {
        if (diagonal == null || superior == null || superior.length != diagonal.length - 1)
            throw new DimensionMatrizException("Bandas incompatibles");
        BigComplex[][] d = datosDimensiones(diagonal.length, diagonal.length, BigComplex.ZERO, false);
        for (int i = 0; i < diagonal.length; i++) {
            d[i][i] = Objects.requireNonNull(diagonal[i]);
            if (i + 1 < diagonal.length) d[i][i + 1] = Objects.requireNonNull(superior[i]);
        }
        return new Matriz(d);
    }
    public static Matriz bidiagonalInferior(BigComplex[] diagonal, BigComplex[] inferior) {
        return bidiagonalSuperior(diagonal, inferior).traspuesta();
    }
    public static Matriz tridiagonal(BigComplex[] inferior, BigComplex[] diagonal, BigComplex[] superior) {
        if (diagonal == null || inferior == null || superior == null
                || inferior.length != diagonal.length - 1 || superior.length != diagonal.length - 1)
            throw new DimensionMatrizException("Bandas incompatibles");
        Matriz m = bidiagonalSuperior(diagonal, superior);
        for (int i = 1; i < diagonal.length; i++) m = m.conElemento(i, i - 1, inferior[i - 1]);
        return m;
    }
    public static Matriz banda(int filas, int columnas, int inferior, int superior, BigComplex valor) {
        if (inferior < 0 || superior < 0) throw new DimensionMatrizException("Ancho de banda negativo");
        BigComplex[][] d = datosDimensiones(filas, columnas, BigComplex.ZERO, false);
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            if (i - j <= inferior && j - i <= superior) d[i][j] = Objects.requireNonNull(valor);
        return new Matriz(d);
    }
    public static Matriz hessenbergSuperior(BigComplex[][] datos) {
        Matriz m = new Matriz(datos);
        if (!m.esHessenbergSuperior()) throw new DimensionMatrizException("No es Hessenberg superior");
        return m;
    }
    public static Matriz desdeMatrizAdyacencia(BigComplex[][] datos) {
        Matriz m = new Matriz(datos);
        if (!m.esCuadrada() || !m.esNoNegativa()) throw new DimensionMatrizException("Adyacencia inválida");
        return m;
    }
    public static Matriz matrizDeIncidencia(int vertices, int[][] aristas, boolean dirigida) {
        if (vertices <= 0 || aristas == null) throw new DimensionMatrizException("Grafo inválido");
        BigComplex[][] d = datosDimensiones(vertices, Math.max(1, aristas.length), BigComplex.ZERO, false);
        for (int e = 0; e < aristas.length; e++) {
            if (aristas[e] == null || aristas[e].length != 2 || aristas[e][0] < 0 || aristas[e][0] >= vertices
                    || aristas[e][1] < 0 || aristas[e][1] >= vertices) throw new DimensionMatrizException("Arista inválida");
            d[aristas[e][0]][e] = dirigida ? Complejos.MINUS_ONE : BigComplex.ONE;
            d[aristas[e][1]][e] = BigComplex.ONE;
        }
        return new Matriz(d);
    }
    public static Matriz matrizCompanera(BigComplex... coeficientesMonicos) {
        if (coeficientesMonicos == null || coeficientesMonicos.length < 2
                || !Complejos.aproximadamenteIguales(coeficientesMonicos[0], BigComplex.ONE,
                Complejos.TOLERANCIA_DEFECTO, Complejos.CONTEXTO_DEFECTO))
            throw new DimensionMatrizException("Se requieren coeficientes mónicos en orden descendente");
        int n = coeficientesMonicos.length - 1;
        BigComplex[][] d = datosDimensiones(n, n, BigComplex.ZERO, false);
        for (int i = 1; i < n; i++) d[i][i - 1] = BigComplex.ONE;
        for (int i = 0; i < n; i++) d[i][n - 1] = coeficientesMonicos[n - i].negate();
        return new Matriz(d);
    }
    public static Matriz bloqueJordan(BigComplex lambda, int n) {
        BigComplex[][] d = datosDimensiones(n, n, BigComplex.ZERO, false);
        for (int i = 0; i < n; i++) {
            d[i][i] = Objects.requireNonNull(lambda);
            if (i + 1 < n) d[i][i + 1] = BigComplex.ONE;
        }
        return new Matriz(d);
    }

    private static BigComplex[][] datosDimensiones(int f, int c, BigComplex valor, boolean diagonal) {
        if (f <= 0 || c <= 0) throw new DimensionMatrizException("Dimensiones positivas requeridas");
        Objects.requireNonNull(valor);
        BigComplex[][] d = new BigComplex[f][c];
        for (int i = 0; i < f; i++) for (int j = 0; j < c; j++) d[i][j] = diagonal ? (i == j ? valor : BigComplex.ZERO) : valor;
        return d;
    }
    private static BigComplex[][] validarYCopiar(BigComplex[][] a, MathContext mc) {
        Objects.requireNonNull(a, "El array no puede ser null");
        if (a.length == 0 || a[0] == null || a[0].length == 0) throw new DimensionMatrizException("Matriz vacía");
        int c = a[0].length; BigComplex[][] d = new BigComplex[a.length][c];
        for (int i = 0; i < a.length; i++) {
            if (a[i] == null || a[i].length != c) throw new DimensionMatrizException("Array irregular");
            for (int j = 0; j < c; j++) d[i][j] = Objects.requireNonNull(a[i][j], "Entrada null").round(mc);
        }
        return d;
    }
    private static BigComplex[][] copiarEstructura(BigComplex[][] a) {
        BigComplex[][] d = new BigComplex[a.length][];
        for (int i = 0; i < a.length; i++) d[i] = a[i].clone();
        return d;
    }
    private static BigComplex[][] copiar(MatrizCompleja m) {
        Objects.requireNonNull(m); BigComplex[][] d = new BigComplex[m.filas()][m.columnas()];
        for (int i = 0; i < m.filas(); i++) for (int j = 0; j < m.columnas(); j++) d[i][j] = m.get(i, j);
        return d;
    }
    private static BigComplex[][] convertir(double[][] a) {
        Objects.requireNonNull(a); BigComplex[][] d = new BigComplex[a.length][];
        for (int i = 0; i < a.length; i++) {
            if (a[i] == null) throw new DimensionMatrizException("Fila null");
            d[i] = new BigComplex[a[i].length];
            for (int j = 0; j < a[i].length; j++) d[i][j] = Complejos.de(a[i][j]);
        }
        return d;
    }
    private static BigComplex[][] convertir(Double[][] a) {
        Objects.requireNonNull(a); BigComplex[][] d = new BigComplex[a.length][];
        for (int i = 0; i < a.length; i++) {
            if (a[i] == null) throw new DimensionMatrizException("Fila null");
            d[i] = new BigComplex[a[i].length];
            for (int j = 0; j < a[i].length; j++) d[i][j] = Complejos.de(Objects.requireNonNull(a[i][j]));
        }
        return d;
    }
    private static BigComplex[][] convertir(BigDecimal[][] a) {
        Objects.requireNonNull(a); BigComplex[][] d = new BigComplex[a.length][];
        for (int i = 0; i < a.length; i++) {
            if (a[i] == null) throw new DimensionMatrizException("Fila null");
            d[i] = new BigComplex[a[i].length];
            for (int j = 0; j < a[i].length; j++) d[i][j] = BigComplex.valueOf(Objects.requireNonNull(a[i][j]));
        }
        return d;
    }
    private static BigComplex[][] convertir(String[][] a) {
        Objects.requireNonNull(a); BigComplex[][] d = new BigComplex[a.length][];
        for (int i = 0; i < a.length; i++) {
            if (a[i] == null) throw new DimensionMatrizException("Fila null");
            d[i] = new BigComplex[a[i].length];
            for (int j = 0; j < a[i].length; j++) d[i][j] = Complejos.de(a[i][j]);
        }
        return d;
    }
}
