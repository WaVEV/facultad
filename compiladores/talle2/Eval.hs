-- | Módulo para evaluar comandos del lenguaje imperativo simple.
module Eval where

import Lang

-- El módulo 'State' exporta las funciones 'get :: State -> String ->
-- Int' y 'set :: State -> String -> Int -> State'; no pueden usar la
-- representación de 'State'.
import State

-- | La denotación de nuestros comandos será un tipo de datos que
-- describa los posibles comportamientos que hemos analizado en el
-- teórico. La variable del input la usamos para poder decir para
-- que variable se está pidiendo input.
data Omega = Term State
           | Abort State
           | Out Int Omega
           | In Iden (Int -> Omega)

-- Notemos que los elementos de Ω no tienen ninguna relación con la
-- entrada o salida; lo que haremos en cambio será luego “interpretar”
-- esas descripciones de comportamientos como computaciones en la 
-- mónada de IO.

-- | Evaluación de Pares
evalPar :: ParExp -> State -> Bool
evalPar e st = case e of
                    (Par e1 e2) -> Tuple (evalInt e1) (evalInt e2)
                    (PVar v) -> get st v

-- | La evaluación de expresiones enteras serán enteros.
evalInt :: IntExpr -> State -> Int
evalInt e st = case e of 
                (Const a) -> a
                (Var a) -> get st a
                (Neg a) -> - evalInt a st
                (IBin op a b) -> f (evalInt a st) (evalInt b st)
                    where f = case op of
                                    Plus -> (\x y -> x + y)
                                    Minus -> (\x y -> x - y)
                                    Times -> (\x y -> x * y)
                                    Div -> (\x y -> x `div` y)
                                    Mod -> (\x y -> x `mod` y)
                (Pri pe) -> e1
                    where Tuple e1 e2 = evalPar pe
                (Seg pe) -> e2
                    where Tuple e1 e2 = evalPar pe


-- | Evaluación de expresiones booleanas en bool.
evalBool :: Assert -> State -> Bool
evalBool e st = case e of 
                    CTrue -> True
                    CFalse -> False
                    (Not a) -> not (evalBool a st)
                    (ABin op a1 a2) -> f (evalBool a1 st) (evalBool a2 st)
                            where f = case op of
                                        And -> \ x y -> x && y
                                        Or -> \ x y -> x || y
                    (IntAss op e0 e1) -> f (evalInt e0 st) (evalInt e1 st)
                            where f = case op of
                                        Eq -> \ x y -> x == y
                                        NEq -> \ x y -> x /= y
                                        Lt -> \ x y -> x < y


-- | La evaluación de un comando es un elemento de Ω; el único
-- elemento que no tiene una representación explícita es ⊥. ¿Por qué?
-- Para definir esta función les pueden ser muy útil las funciones que
-- están debajo.
eval :: Comm -> State -> Omega
eval (While b c) st = fix fwhile st
    where fwhile w st' | evalBool b st' = (star w) (eval c st')
                       | otherwise =  Term st'
eval c st = case c of
                Skip -> Term st
                (Assign var e) -> Term (set st var value) 
                    where value = evalInt e st
                (If b c0 c1) | evalBool b st -> eval c0 st
                             | otherwise -> eval c1 st
                (Seq c0 c1) -> star (eval c1) (eval c0 st)
                (Newvar var e c) -> (dagger (\o -> set o var (get st var))) (eval c st')
                    where st' = set st var (evalInt e st)
                (Newvar var e1 e2 c) -> (dagger (\o -> set o var(get st var))) (eval c st')
                    where st' = set st var (Tuple (evalInt e1 st) (evalInt e2 st))
                Morite -> Abort st
                (Agarrame c0 c1) -> pluss (eval c1) (eval c0 st)

                (Toma (Par e1 e2)) -> Out n1 (Out n2 (Term st))
                    where Tuple n1 n2 = evalPar (Par e1 e2)
                (Toma (PVar var)) -> Out n1 (Out n2 (Term st))
                    where Tuple n1 n2 = evalPar (Par e1 e2)

                (Toma e) -> Out n (Term st)
                    where n = evalInt e st

                (Dame var) -> In var (\n -> Term (set st var n) )
                (DamePar var) -> In var (\ n1 n2 -> Term(set st var (Tuple n1 n2)))
                (Assign var e1 e2) -> Term(set st var (Tuple e1 e2))



pluss :: (State -> Omega) -> Omega -> Omega
pluss f (Term st) = Term st
pluss f (Abort st) = f st
pluss f (Out n w ) = Out n (pluss f w)
pluss f (In v g) = In v (pluss f . g)

-- | La función star (★) extiende funciones de Σ → Ω a funciones de Ω → Ω. 
-- ¿Por qué no propagamos ⊥? ¿O sí lo estamos propagando?
star :: (State -> Omega) -> Omega -> Omega
star f (Term st) = f st
star f (Abort st) = Abort st
star f (Out n w) = Out n (star f w)
star f (In v g) = In v (star f . g) -- esto es equivalente a \n -> star f (g n)

-- | La función dagger (†) extiende funciones de Σ → Σ a funciones de
-- Ω → Ω, aplicando la función al estado resultante en caso de terminación
-- normal o anormal y en los casos de input y output se comporta similar a
-- star.
dagger :: (State -> State) -> Omega -> Omega
dagger f (Term st) = Term (f st)
dagger f (Abort st) = Abort (f st)
dagger f (Out n w) = Out n (dagger f w)
dagger f (In v g) = In v (dagger f . g) 

-- | En Haskell podemos definir una función que resembla el operador Y del
-- menor punto fijo.
fix :: (d -> d) -> d
fix f = f (fix f)

-- | La siguiente función de alto-orden tiene como menor punto fijo el
--  factorial, para los positivos...
fact :: (Int -> Int) -> (Int -> Int)
fact f 0 = 1
fact f n = n * fact f (n-1)

-- | y de hecho podemos definirlo como tal.
factorial :: Int -> Int
factorial = fix fact

{- 
Para convencerse de que funciona podemos pensar cómo se evalúa la
expresión @factorial 4@:

factorial 4 = (fix fact) 4
→ 
fact (fix fact) 4
→
4 * fact (fix fact) 3
→
4 * (3 * fact (fix fact) 2)
→ 
4 * (3 * (2 * fact (fix fact) 1))
→
4 * (3 * (2 * (1 * fact (fix fact) 0)))
→
4 * (3 * (2 * (1 * 1)))
→→
24
-}
