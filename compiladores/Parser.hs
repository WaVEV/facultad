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
           else fail "ola k ase"

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

pCatch :: TParser Comm
pCatch = do _ <- silentKw KCatchin
            c1 <- pComm
            _ <- silentKw In
            _ <- silentSym OpenCurly 
            c2 <- pComm
            _ <- silentSym CloseCurly
            return $ Catchin c1 c2

pFail = do _ <- silentKw KFail
           return $ Fail

pOutput = do _ <- silentSym Exclamation
             e <- pIntExp
             return $ Output e

pInput = do _ <- silentSym Question
            v <- pvar
            return $ Input v


-- | Parser para comandos que no tienen secuencias. (COMPLETAR)
pComm' :: TParser Comm
pComm' = choice [pAssign,pif,pskip,pNewVar,pwhile, pCatch, pFail,pOutput,pInput]


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
pIntExp = try pIntExpRSumSub <|> try pIntExpRTMD <|> try pBasicIntExp


pIntExpRSumSub:: TParser IntExpr
pIntExpRSumSub = try parse_sum_sub   <|> 
                 try pIntExpRTMD     <|> 
                 try pBasicIntExp    <|>
                 error "Expresión inválida IntExp "
              
    where parse_sum_sub = do    n <- pIntExpRTMD
                                op <- ptk isSumSubOp
                                m <- pIntExp
                                return $ Bexpr n op m
          isSumSubOp (Op TokPlus) = Just Plus
          isSumSubOp (Op TokMinus) = Just Minus
          isSumSubOp _ = Nothing

pIntExpRTMD:: TParser IntExpr
pIntExpRTMD = try parse_times_div <|>
              try pBasicIntExp    
          
    where parse_times_div = do  n <- pBasicIntExp
                                op <- ptk isTimesDivOp
                                m <- pIntExpRTMD
                                return $ Bexpr n op m
          isTimesDivOp (Op TokTimes) = Just Times
          isTimesDivOp (Op TokMod) = Just Mod
          isTimesDivOp (Op TokDiv) = Just Div
          isTimesDivOp _ = Nothing

pBasicIntExp:: TParser IntExpr
pBasicIntExp = try var'         <|>
               try num          <|>
               try not          <|>
               try nest_IntExpr <|>
               error "Expresión inválida IntExp 2"
          
    where var' = do             v <- pvar
                                return $ Var v
          num = do              n <- ptk isNum
                                return $ Const n
          isNum (Num n) = Just n
          isNum _ = Nothing

          nest_IntExpr = do     _ <- silentSym OpenParen
                                n <-  pIntExp
                                _ <- silentSym CloseParen
                                return $ n
          not = do _  <- ptk isNeg
                   n <- pBasicIntExp
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
pAssert = try pBoolConst   <|> 
          try pAssertIFF   <|>
          try pAssertAO    <|>
          try pAssertBasic 

pAssertIFF :: TParser Assert
pAssertIFF = try parse_Iff <|>
             try parse_If 

        where parse_Iff = do    a_1 <- pAssertAO
                                r <- ptk isIff
                                a_2 <- pAssert
                                return $ ABin a_1 r a_2
              isIff (Op TokIff) = Just Iff
              isIff _ = Nothing

              parse_If = do     a_1 <-  pAssertAO
                                r <- ptk isIf
                                a_2 <- pAssert
                                return $ ABin a_1 r a_2
              isIf (Op TokImpl) = Just Impl
              isIf _ = Nothing

pAssertAO :: TParser Assert
pAssertAO = try parse_Or  <|>
            try parse_And <|>
            try pAssertBasic

        where parse_Or = do     a_1 <-  pAssertBasic
                                r <- ptk isOr
                                a_2 <- pAssertAO
                                return $ ABin a_1 r a_2
              isOr (Op TokOr) = Just Or
              isOr _ = Nothing

              parse_And = do    a_1 <-  pAssertBasic
                                r <- ptk isAnd
                                a_2 <- pAssertAO
                                return $ ABin a_1 r a_2
              isAnd (Op TokAnd) = Just And
              isAnd _ = Nothing

pAssertBasic :: TParser Assert
pAssertBasic = try pBoolConst <|>
               try nest_assert <|>
               try rel_assert  <|>
               try not         <|>
               try rel_assert 


        where nest_assert = do  _ <- silentSym OpenParen 
                                st <- pAssert
                                _ <- silentSym CloseParen
                                return $ st

              not = do          _ <- ptk isNeg
                                st <- pAssertBasic
                                return $ Not st
              isNeg (Op TokNot) = Just Not
              isNeg _ = Nothing

              rel_assert = do   n <- pIntExp
                                r <- ptk isOpRel
                                m <- pIntExp
                                return $ Relacion n r m

              isOpRel (Op TokEqual) = Just Eq
              isOpRel (Op TokNotEq) = Just NEq
              isOpRel (Op TokLess) = Just Less
              isOpRel (Op TokLtEq) = Just LtEq
              isOpRel (Op TokGreater) = Just Greater
              isOpRel (Op TokGtEq) = Just GtEq
              isOpRel _ = Nothing


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

