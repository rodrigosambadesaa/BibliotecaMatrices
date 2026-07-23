package comun.resultados;

import comun.MatrizCompleja;

import java.math.BigDecimal;
import java.util.List;

/** Resultado inmutable de un sistema, incluyendo clasificación y residuo. */
public record ResultadoSistema<M extends MatrizCompleja>(
        TipoSolucion tipo, M solucionParticular, List<M> baseHomogenea, BigDecimal residuo) {
    public ResultadoSistema {
        baseHomogenea = List.copyOf(baseHomogenea);
    }
}
