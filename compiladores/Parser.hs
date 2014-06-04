-- | Parser para el lenguaje imperativo simple. El parser
-- toma como entrada una lista de tokens.

module Parser where

import qualified Control.Exception as E
import Lang
import Lexer


import Text.Parsec

type TParser = Parsec [Token] ()

parser :: [Token] -> Either ParseError Comm
parser = runParser pcomando () "lexer"

pvar' = pvar >>= \v -> empty >> return v

pcomando = pComm >>= \v -> empty >> return v


empty :: TParser ()
empty = do s <- getInput
           if null s then return ()
           else fail "FAIL EN EMPTY"

-- | Un parser muy sencillo para skip
pskip :: TParser Comm
pskip = do _ <- silentKw KSkip
           return Skip

-- | Un parser para if ... then ... else 
pif :: TParser Comm
pif = do _ <- silentKw KIf
         b <- pAssert
         _ <- silentKw Then
         c1 <- pComm
         _ <- silentKw Else
         c2 <- pComm
         _ <- silentKw Fi
         return $ If b c1 c2

pwhile :: TParser Comm
pwhile = do _ <- silentKw KWhile
            b <- pAssert
            _ <- silentKw Do
            c <- pComm
            _ <- silentKw Od
            return $ While b c

pAssign :: TParser Comm
pAssign = do v <- pvar
             _ <- silentOp TokAssign
             e <- pIntExp
             return $ Assign v e


pNewVar :: TParser Comm
pNewVar = do _ <- silentKw KNewVar
             v <- pvar
             _ <- silentOp TokAssign
             e <- pIntExp
             _ <- silentKw In
             _ <- silentSym OpenCurly 
             c <- pComm
             _ <- silentSym CloseCurly
             return $ NewVar v e c

pOutput :: TParser Comm
pOutput = do _ <- silentSym Admiration
             v <- pIntExp
             return $ Output v

pFail :: TParser Comm
pFail = do _ <- silentKw KFail
           return $ Fail


pInput:: TParser Comm
pInput = do _ <- silentSym Interrogation
            v <- pvar
            return $ Input v


pCatch::TParser Comm
pCatch = do _ <- silentKw Catchin
            v <- pComm
            _ <- silentKw In
            _ <- silentSym OpenCurly 
            c <- pComm
            _ <- silentSym CloseCurly
            return $ Catch v c



-- | Parser para comandos que no tienen secuencias. (COMPLETAR)
pComm' :: TParser Comm
pComm' = choice [pAssign,pif,pskip,pNewVar,pwhile,pFail,pOutput,pInput,pCatch]


-- | Parser de secuencias
pseq :: TParser Comm
pseq = do c <- pComm'
          _ <- silentOp TokSemicolon
          c' <- pComm
          return $ Seq c c'

-- | Parser de comandos: o bien una secuencia o bien un comando solito.
pComm :: TParser Comm
pComm = try pseq <|> pComm'


pvar :: TParser Iden
pvar = ptk str
    where str (Str v) = Just v
          str _ = Nothing




pIntExp:: TParser IntExpr
pIntExp = ie_plus 
              
{-
ie_plus = ie_times ("+","-") ie_plus | ie_times
ie_times = ie_0 ("*","%","/") ie_times | ie_0
ie_0 = var | -ie_plus | const | (ie_plus)
-}

ie_plus:: TParser IntExpr
ie_plus = try parse_sum <|> try ie_times <|> error "Expresión inválida 1 "
    where parse_sum = do
                      n <-  ie_times
                      op <- ptk esOperador
                      m <- ie_plus
                      return $ Bexpr n op m
          esOperador (Op TokPlus) = Just Plus
          esOperador (Op TokMinus) = Just Minus
          esOperador _ = Nothing


ie_times:: TParser IntExpr
ie_times = try parse_times <|> try  ie_zero <|> error "Expresión inválida 2"
    where parse_times = do
                        n <-  ie_zero
                        op <- ptk esOperador
                        m <- ie_times
                        return $ Bexpr n op m
          esOperador (Op TokTimes) = Just Times
          esOperador (Op TokMod) = Just Mod
          esOperador (Op TokDiv) = Just Div
          esOperador _ = Nothing






ie_zero:: TParser IntExpr
ie_zero = try var' <|>  try negacion_general <|> try num <|> try  parentesis <|> error "Expresión inválida 3"
            where var' = do 
                            v <- pvar
                            return $ Var v
                  num = do n <- ptk isNum
                           return $ Const n
                  isNum (Num n) = Just n
                  isNum _ = Nothing
                  parentesis = do 
                                    _ <- silentSym OpenParen
                                    n <-  pIntExp
                                    _ <- silentSym CloseParen
                                    return $ n
                  negacion_general = do _  <- ptk isNeg
                                        n <- pIntExp
                                        return $ Neg n
                  isNeg (Op TokMinus) = Just Neg
                  isNeg _ = Nothing
                



-- | Parser para constantes booleanas.
pBoolConst :: TParser Assert
pBoolConst = do b <- pkws [ (KTrue,CTrue) 
                          , (KFalse,CFalse)
                          ]
                return b

 

pAssert :: TParser Assert
pAssert = p_Sii

p_Sii :: TParser Assert
p_Sii = try parse_Sii <|> try  p_Impl <|>  error "p_sii assert"
        where  
              parse_Sii = do       a_1 <-  p_Impl
                                   r <- ptk esIff
                                   a_2 <- p_Sii
                                   return $ ABin a_1 r a_2
              esIff (Op TokIff) = Just Iff
              esIff _ = Nothing


p_Impl :: TParser Assert
p_Impl = try parse_Impl <|> try  p_Or <|>  error "p_impl assert"
         where  
               parse_Impl = do      a_1 <-  p_Or
                                    r <- ptk esMultiple
                                    a_2 <- p_Impl
                                    return $ ABin a_1 r a_2
               esMultiple (Op TokImpl) = Just Impl
               esMultiple _ = Nothing


p_Or :: TParser Assert
p_Or = try parse_Or <|> try p_And <|>  error "p_or assert"
       where  
             parse_Or = do      a_1 <-  p_And
                                r <- ptk esMultiple
                                a_2 <- p_Or
                                return $ ABin a_1 r a_2
             esMultiple (Op TokOr) = Just Or
             esMultiple _ = Nothing


p_And :: TParser Assert
p_And = try parse_And <|> try p_zero <|>  error "p_And assert"
        where  
              parse_And = do      a_1 <-  p_zero
                                  r <- ptk esMultiple
                                  a_2 <- p_And
                                  return $ ABin a_1 r a_2
              esMultiple (Op TokAnd) = Just And
              esMultiple _ = Nothing

p_zero :: TParser Assert
p_zero = try pBoolConst <|> try negacion <|> try parentesis_ass  <|> try relaciones <|> error "p_zero assert"
         where           
               parentesis_ass = do _ <- silentSym OpenParen 
                                   st <- pAssert
                                   _ <- silentSym CloseParen
                                   return $ st

               negacion = do  _ <- ptk esNegacion
                              st <- pAssert
                              return $ Not st
               esNegacion (Op TokNot) = Just Not
               esNegacion _ = Nothing

               relaciones = do 
                                n <-  pIntExp
                                r <- ptk esOperadorRelacional
                                m <- pIntExp
                                return $ Relacion n r m

               esOperadorRelacional (Op TokEqual) = Just Eq
               esOperadorRelacional (Op TokNotEq) = Just NEq
               esOperadorRelacional (Op TokLess) = Just Less
               esOperadorRelacional (Op TokLtEq) = Just LtEq
               esOperadorRelacional (Op TokGreater) = Just Greater
               esOperadorRelacional (Op TokGtEq) = Just GtEq
               esOperadorRelacional _ = Nothing


-- | Consumimos un keyword y nada más.
silentKw :: Keyword -> TParser ()
silentKw k = silent isK k
    where isK (Kw k') = Just k'
          isK _ = Nothing

-- | Consumimos un operador y nada más.
silentOp :: Operator -> TParser ()
silentOp o = silent isOp o
    where isOp (Op o') = Just o'
          isOp _ = Nothing

-- | Consumimos un símbolo y nada más.
silentSym :: Symbol -> TParser ()
silentSym s = silent isSym s
    where isSym (Sym s') = Just s'
          isSym _ = Nothing

-- | Si tenemos la suerte que isA para el token actual devuelva justo
-- a', entonces comparamos a con a', en ese caso consumimos el token
-- actual si no, no consumimos nada.
silent :: Eq a => (Token -> Maybe a) -> a -> TParser ()
silent isA a = ptk check
    where check t = case isA t of
                      Just a' -> if a == a' then Just () else Nothing
                      Nothing -> Nothing

-- | 
pkws :: [(Keyword,a)] -> TParser a
pkws cs = ptk l
    where l (Kw k) = lookup k cs
          l _ = Nothing

-- | Dada una función que, con suerte, nos devuelve un a, vemos 
-- un token, y si tuvimos suerte lo consumimos y devolvemos un a.
ptk :: (Token -> Maybe a) -> TParser a
ptk f = tokenPrim show incpos f

-- | Cada vez que consumimos un token incrementamos la posición de
-- nuestro parser (notar que ya perdimos la relación con la posición
-- en el string original).
incpos :: SourcePos -> Token -> [Token] -> SourcePos
incpos p _ _ = incSourceColumn p 1

