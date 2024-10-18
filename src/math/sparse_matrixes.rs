use std::ops::Index;
use std::ops::IndexMut;

struct TridiagonalSparseMatrix<NumberType> {
    subdiagonal: Vec<NumberType>,
    maindiagonal: Vec<NumberType>,
    updiagonal: Vec<NumberType>,
}

impl<NumberType: std::convert::From<i32> + Clone> TridiagonalSparseMatrix<NumberType> {
    fn new(n: usize) -> TridiagonalSparseMatrix<NumberType> {
        TridiagonalSparseMatrix {
            subdiagonal: vec![NumberType::from(0); n],
            maindiagonal: vec![NumberType::from(0); n],
            updiagonal: vec![NumberType::from(0); n],
        }
    }
}

impl<'a, NumberType> Index<(usize, usize)> for TridiagonalSparseMatrix<NumberType> {
    type Output = NumberType;
    fn index(&self, idx: (usize, usize)) -> &NumberType {
        let (x, y) = idx;
        match (x, y) {
            (x, y) if x == y => self.maindiagonal.get(x).unwrap(),
            (x, y) if x == y - 1 => self.subdiagonal.get(x).unwrap(),
            (x, y) if x == y + 1 => self.updiagonal.get(x).unwrap(),
            _ => panic!("Index out of bounds"),
        }
    }
}

impl<'a, NumberType> IndexMut<(usize, usize)> for TridiagonalSparseMatrix<NumberType> {
    fn index_mut(&mut self, idx: (usize, usize)) -> &mut NumberType {
        let (x, y) = idx;
        match (x, y) {
            (x, y) if x == y => self.maindiagonal.get_mut(x).unwrap(),
            (x, y) if x == y - 1 => self.subdiagonal.get_mut(x).unwrap(),
            (x, y) if x == y + 1 => self.updiagonal.get_mut(x).unwrap(),
            _ => panic!("Index out of bounds"),
        }
    }
}
