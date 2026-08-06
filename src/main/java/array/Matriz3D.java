package array;

import ch.obermuhlner.math.big.BigComplex;
import comun.AbstractMatriz3D;
import comun.Complejos;
import comun.Matriz3DCompleja;
import comun.MatrizCompleja;
import comun.excepciones.DimensionMatrizException;
import comun.excepciones.OperacionNoDefinidaException;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Objects;
import java.util.Random;

/** Matriz 3D inmutable almacenada realmente en {@code BigComplex[][][]}. */
public final class Matriz3D extends AbstractMatriz3D<Matriz3D, Matriz> {
    private final BigComplex[][][] matriz;

    public Matriz3D(int capas, int filas, int columnas) {
        this(capas, filas, columnas, BigComplex.ZERO);
    }

    public Matriz3D(int capas, int filas, int columnas, BigComplex valorInicial) {
        this(nuevosDatos(capas, filas, columnas, Objects.requireNonNull(valorInicial)),
                Complejos.CONTEXTO_DEFECTO, Complejos.TOLERANCIA_DEFECTO);
    }

    public Matriz3D(int capas, int filas, int columnas, MathContext contexto) {
        this(nuevosDatos(capas, filas, columnas, BigComplex.ZERO), contexto,
                BigDecimal.ONE.scaleByPowerOfTen(-Math.max(1, contexto.getPrecision() - 6)));
    }

    public Matriz3D(double[][][] datos) { this(convertir(datos)); }
    public Matriz3D(Double[][][] datos) { this(convertir(datos)); }
    public Matriz3D(BigDecimal[][][] datos) { this(convertir(datos)); }
    public Matriz3D(String[][][] datos) { this(convertir(datos)); }
    public Matriz3D(BigComplex[][][] datos) {
        this(datos, Complejos.CONTEXTO_DEFECTO, Complejos.TOLERANCIA_DEFECTO);
    }

    public Matriz3D(Matriz3D otra) {
        this((Matriz3DCompleja) Objects.requireNonNull(otra));
    }

    public Matriz3D(Matriz3DCompleja otra) {
        this(otra, otra.mathContext(), otra.tolerancia());
    }

    public Matriz3D(Matriz3DCompleja otra, MathContext contexto, BigDecimal tolerancia) {
        this(copiar(otra), contexto, tolerancia);
    }

    public Matriz3D(BigComplex[][][] datos, MathContext contexto, BigDecimal tolerancia) {
        this(validarYCopiar(datos, contexto), contexto, tolerancia, true);
    }

    private Matriz3D(BigComplex[][][] datos, MathContext contexto, BigDecimal tolerancia, boolean validados) {
        super(datos.length, datos[0].length, datos[0][0].length, contexto, tolerancia);
        matriz = validados ? copiarEstructura(datos) : validarYCopiar(datos, contexto);
    }

    @Override public BigComplex get(int capa, int fila, int columna) {
        validarIndice(capa, fila, columna);
        return matriz[capa][fila][columna];
    }

    @Override protected Matriz3D crear(BigComplex[][][] datos, MathContext contexto, BigDecimal tolerancia) {
        return new Matriz3D(datos, contexto, tolerancia);
    }

    @Override protected Matriz crearPlano(BigComplex[][] datos, MathContext contexto, BigDecimal tolerancia) {
        return new Matriz(datos, contexto, tolerancia);
    }

    public static Matriz3D cero(int capas, int filas, int columnas) { return new Matriz3D(capas, filas, columnas); }

    public static Matriz3D rellena(int capas, int filas, int columnas, BigComplex valor) {
        return new Matriz3D(capas, filas, columnas, valor);
    }

    /** Crea varias matrices identidad rectangulares, una en cada capa. */
    public static Matriz3D identidadPorCapas(int capas, int filas, int columnas) {
        BigComplex[][][] datos = nuevosDatos(capas, filas, columnas, BigComplex.ZERO);
        for (int k = 0; k < capas; k++)
            for (int i = 0; i < Math.min(filas, columnas); i++) datos[k][i][i] = BigComplex.ONE;
        return new Matriz3D(datos);
    }

    public static Matriz3D desdeMatriz(MatrizCompleja matriz) {
        return desdeCapas(matriz);
    }

    public static Matriz3D desdeCapas(MatrizCompleja... capas) {
        Objects.requireNonNull(capas, "Las capas no pueden ser null");
        if (capas.length == 0 || capas[0] == null) throw new DimensionMatrizException("Se requiere alguna capa");
        MatrizCompleja primera = capas[0];
        BigComplex[][][] datos = new BigComplex[capas.length][primera.filas()][primera.columnas()];
        for (int k = 0; k < capas.length; k++) {
            MatrizCompleja capa = Objects.requireNonNull(capas[k], "Una capa no puede ser null");
            if (capa.filas() != primera.filas() || capa.columnas() != primera.columnas())
                throw new DimensionMatrizException("Las capas deben tener las mismas dimensiones");
            if (!capa.mathContext().equals(primera.mathContext()) || capa.tolerancia().compareTo(primera.tolerancia()) != 0)
                throw new OperacionNoDefinidaException("Las capas deben compartir política de precisión");
            for (int i = 0; i < primera.filas(); i++)
                for (int j = 0; j < primera.columnas(); j++) datos[k][i][j] = capa.get(i, j);
        }
        return new Matriz3D(datos, primera.mathContext(), primera.tolerancia());
    }

    public static Matriz3D aleatoria(int capas, int filas, int columnas, long semilla) {
        Random random = new Random(semilla);
        BigComplex[][][] datos = nuevosDatos(capas, filas, columnas, BigComplex.ZERO);
        for (int k = 0; k < capas; k++) for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            datos[k][i][j] = Complejos.de(random.nextDouble() * 2 - 1);
        return new Matriz3D(datos);
    }

    private static BigComplex[][][] validarYCopiar(BigComplex[][][] origen, MathContext contexto) {
        Objects.requireNonNull(origen, "Los datos no pueden ser null");
        Objects.requireNonNull(contexto, "El MathContext no puede ser null");
        if (origen.length == 0 || origen[0] == null || origen[0].length == 0
                || origen[0][0] == null || origen[0][0].length == 0)
            throw new DimensionMatrizException("Matriz 3D vacía");
        int filas = origen[0].length;
        int columnas = origen[0][0].length;
        BigComplex[][][] copia = new BigComplex[origen.length][filas][columnas];
        for (int k = 0; k < origen.length; k++) {
            if (origen[k] == null || origen[k].length != filas) throw new DimensionMatrizException("Capas irregulares");
            for (int i = 0; i < filas; i++) {
                if (origen[k][i] == null || origen[k][i].length != columnas)
                    throw new DimensionMatrizException("Filas irregulares en la capa " + k);
                for (int j = 0; j < columnas; j++)
                    copia[k][i][j] = Objects.requireNonNull(origen[k][i][j], "Un elemento no puede ser null").round(contexto);
            }
        }
        return copia;
    }

    private static BigComplex[][][] copiarEstructura(BigComplex[][][] origen) {
        BigComplex[][][] copia = new BigComplex[origen.length][][];
        for (int k = 0; k < origen.length; k++) {
            copia[k] = new BigComplex[origen[k].length][];
            for (int i = 0; i < origen[k].length; i++) copia[k][i] = origen[k][i].clone();
        }
        return copia;
    }

    private static BigComplex[][][] copiar(Matriz3DCompleja origen) {
        Objects.requireNonNull(origen, "La matriz 3D no puede ser null");
        BigComplex[][][] datos = new BigComplex[origen.capas()][origen.filas()][origen.columnas()];
        for (int k = 0; k < origen.capas(); k++) for (int i = 0; i < origen.filas(); i++)
            for (int j = 0; j < origen.columnas(); j++) datos[k][i][j] = origen.get(k, i, j);
        return datos;
    }

    private static BigComplex[][][] convertir(double[][][] origen) {
        Objects.requireNonNull(origen); BigComplex[][][] datos = new BigComplex[origen.length][][];
        for (int k = 0; k < origen.length; k++) {
            if (origen[k] == null) throw new DimensionMatrizException("Capa null");
            datos[k] = new BigComplex[origen[k].length][];
            for (int i = 0; i < origen[k].length; i++) {
                if (origen[k][i] == null) throw new DimensionMatrizException("Fila null");
                datos[k][i] = new BigComplex[origen[k][i].length];
                for (int j = 0; j < origen[k][i].length; j++) datos[k][i][j] = Complejos.de(origen[k][i][j]);
            }
        }
        return datos;
    }

    private static BigComplex[][][] convertir(Double[][][] origen) {
        Objects.requireNonNull(origen); BigComplex[][][] datos = new BigComplex[origen.length][][];
        for (int k = 0; k < origen.length; k++) {
            if (origen[k] == null) throw new DimensionMatrizException("Capa null");
            datos[k] = new BigComplex[origen[k].length][];
            for (int i = 0; i < origen[k].length; i++) {
                if (origen[k][i] == null) throw new DimensionMatrizException("Fila null");
                datos[k][i] = new BigComplex[origen[k][i].length];
                for (int j = 0; j < origen[k][i].length; j++) datos[k][i][j] = Complejos.de(Objects.requireNonNull(origen[k][i][j]));
            }
        }
        return datos;
    }

    private static BigComplex[][][] convertir(BigDecimal[][][] origen) {
        Objects.requireNonNull(origen); BigComplex[][][] datos = new BigComplex[origen.length][][];
        for (int k = 0; k < origen.length; k++) {
            if (origen[k] == null) throw new DimensionMatrizException("Capa null");
            datos[k] = new BigComplex[origen[k].length][];
            for (int i = 0; i < origen[k].length; i++) {
                if (origen[k][i] == null) throw new DimensionMatrizException("Fila null");
                datos[k][i] = new BigComplex[origen[k][i].length];
                for (int j = 0; j < origen[k][i].length; j++) datos[k][i][j] = BigComplex.valueOf(Objects.requireNonNull(origen[k][i][j]));
            }
        }
        return datos;
    }

    private static BigComplex[][][] convertir(String[][][] origen) {
        Objects.requireNonNull(origen); BigComplex[][][] datos = new BigComplex[origen.length][][];
        for (int k = 0; k < origen.length; k++) {
            if (origen[k] == null) throw new DimensionMatrizException("Capa null");
            datos[k] = new BigComplex[origen[k].length][];
            for (int i = 0; i < origen[k].length; i++) {
                if (origen[k][i] == null) throw new DimensionMatrizException("Fila null");
                datos[k][i] = new BigComplex[origen[k][i].length];
                for (int j = 0; j < origen[k][i].length; j++) datos[k][i][j] = Complejos.de(Objects.requireNonNull(origen[k][i][j]));
            }
        }
        return datos;
    }
}
