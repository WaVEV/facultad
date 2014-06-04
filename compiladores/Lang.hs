-- | Definición del lenguaje.
module Lang where

-- | Los nombres de variables son strings, pero quizás 
-- querramos cambiarlo luego.
type Iden = String

-- | Expresiones enteras COMPLETAR
data IntExpr = Const Int
             | Var Iden
             | Neg IntExpr
             | Bexpr IntExpr OpInt IntExpr
             deriving Show



-- | Operadores binarios enteros, `Div` corresponde al cociente.
data OpInt = Plus | Minus | Times | Div | Mod
           deriving Show

-- | Expresiones booleanas COMPLETAR
data Assert = CTrue
            | CFalse
            | Not Assert
            | ABin Assert OpBool Assert
            | Relacion IntExpr OpRel IntExpr
              deriving Show


-- | Operadores de relación (COMPLETAR)
data OpRel = Eq | NEq |Less | Greater |GtEq | Equal | LtEq
           deriving Show

-- | Operadores binarios booleanos (COMPLETAR)
data OpBool = And | Or | Impl | Iff
            deriving Show

-- | Comandos (COMPLETAR)

data Comm = Skip 
          | Assign Iden IntExpr
          | If Assert Comm Comm
          | Seq Comm Comm
          | While Assert Comm
          | Fail
          | Catchin Comm Comm 
          | Output IntExpr
          | Input String
          | NewVar String IntExpr Comm
            deriving Show



-- Adivinen por qué están comentadas estas cosas.
-- -- | Nombres cómodos (aka combinadores) para construcciones:
-- plus :: IntExpr -> IntExpr -> IntExpr
-- plus = IBin Plus

-- infixr 8 .+
-- -- | Versión con operadores infijos
-- (.+) :: IntExpr -> IntExpr -> IntExpr
-- (.+) = IBin Plus
