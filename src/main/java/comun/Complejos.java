package comun;

import ch.obermuhlner.math.big.BigComplex;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Objects;

/**
 * Operaciones centrales sobre complejos. Todas las conversiones decimales usan
 * texto decimal (también para {@code double}) y todo redondeo es explícito.
 */
public final class Complejos {
    public static final BigComplex ZERO = BigComplex.ZERO;
    public static final BigComplex ONE = BigComplex.ONE;
    public static final BigComplex I = BigComplex.I;
    public static final BigComplex MINUS_ONE = BigComplex.valueOf(BigDecimal.ONE.negate());
    public static final MathContext CONTEXTO_DEFECTO = MathContext.DECIMAL128;
    public static final BigDecimal TOLERANCIA_DEFECTO = new BigDecimal("1E-28");

    private Complejos() {
    }

    public static BigComplex de(Number numero) {
        Objects.requireNonNull(numero, "El número no puede ser null");
        if (numero instanceof BigDecimal decimal) {
            return BigComplex.valueOf(decimal);
        }
        return BigComplex.valueOf(new BigDecimal(numero.toString()));
    }

    public static BigComplex de(String texto) {
        Objects.requireNonNull(texto, "El texto no puede ser null");
        String s = texto.trim().replace(" ", "");
        if (s.endsWith("i")) {
            String base = s.substring(0, s.length() - 1);
            int separador = -1;
            for (int i = 1; i < base.length(); i++) {
                char c = base.charAt(i);
                if ((c == '+' || c == '-') && base.charAt(i - 1) != 'E' && base.charAt(i - 1) != 'e') {
                    separador = i;
                }
            }
            if (separador >= 0) {
                return BigComplex.valueOf(new BigDecimal(base.substring(0, separador)),
                        decimalImaginario(base.substring(separador)));
            }
            return BigComplex.valueOf(BigDecimal.ZERO, decimalImaginario(base));
        }
        return BigComplex.valueOf(new BigDecimal(s));
    }

    private static BigDecimal decimalImaginario(String s) {
        return switch (s) {
            case "", "+" -> BigDecimal.ONE;
            case "-" -> BigDecimal.ONE.negate();
            default -> new BigDecimal(s);
        };
    }

    public static BigComplex redondear(BigComplex z, MathContext mc) {
        return Objects.requireNonNull(z, "El complejo no puede ser null").round(mc);
    }

    public static BigDecimal modulo(BigComplex z, MathContext mc) {
        return Objects.requireNonNull(z, "El complejo no puede ser null").abs(mc);
    }

    public static BigComplex conjugado(BigComplex z) {
        return Objects.requireNonNull(z, "El complejo no puede ser null").conjugate();
    }

    public static boolean esCero(BigComplex z, BigDecimal tolerancia, MathContext mc) {
        validarTolerancia(tolerancia);
        return modulo(z, mc).compareTo(tolerancia) <= 0;
    }

    public static boolean aproximadamenteIguales(BigComplex a, BigComplex b,
                                                  BigDecimal tolerancia, MathContext mc) {
        return esCero(a.subtract(b, mc), tolerancia, mc);
    }

    public static void validarTolerancia(BigDecimal tolerancia) {
        Objects.requireNonNull(tolerancia, "La tolerancia no puede ser null");
        if (tolerancia.signum() < 0) {
            throw new IllegalArgumentException("La tolerancia no puede ser negativa");
        }
    }
}
