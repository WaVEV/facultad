module State (State, DosInt, get,set,initState,iota_int,iota_tuple,primero,segundo,intToDosInt,iota1_TupleInt,iota2_TupleInt) where

newtype DosInt = DosInt (Int, Int)
instance Show DosInt where
	show (DosInt a) = show a

data TupleInt = L1 Int | L2 DosInt

instance Show TupleInt where
	show (L1 a) = (show a)
	show (L2 a) = (show a)

newtype State = State [(String, TupleInt)]

iota_int :: Int -> TupleInt
iota_int a = L1 a

iota_tuple :: DosInt -> TupleInt
iota_tuple a = L2 a

iota1_TupleInt :: TupleInt -> Int
iota1_TupleInt (L1 a) = a
iota1_TupleInt _ = error "Dominio Incorrecto"

iota2_TupleInt :: TupleInt -> DosInt
iota2_TupleInt (L2 a) = a
iota2_TupleInt _ = error "Dominio Incorrecto"

primero:: DosInt -> Int
primero (DosInt a) = fst a

segundo:: DosInt -> Int
segundo (DosInt a) = snd a

myid :: TupleInt -> TupleInt
myid a = a

intToDosInt:: Int->Int->DosInt
intToDosInt a b = DosInt (a, b)

get :: State -> String -> TupleInt
get (State s) v = maybe (iota_int 0) myid $ lookup v s

set :: State -> String -> TupleInt -> State
set (State s) v i = State $ update s v i
  where update [] v i = [(v,i)]
        update ((v',j):vs) v i | v == v' = (v,i):vs
                               | v /= v' = (v',j):update vs v i

instance Show State where
	show (State s) = show (take 4 s) ++ "…"

initState :: State
initState = State []
