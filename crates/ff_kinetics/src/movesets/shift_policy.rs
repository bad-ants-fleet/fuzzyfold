
pub trait ShiftPolicy: Copy {
    #[inline(always)]
    fn three_way() -> bool { false }

    #[inline(always)]
    fn four_way() -> bool { false }
}

#[derive(Debug, Clone, Copy)]
pub struct NoShift;

#[derive(Debug, Clone, Copy)]
pub struct ThreeWayOnly;

#[derive(Debug, Clone, Copy)]
pub struct FourWayOnly;

#[derive(Debug, Clone, Copy)]
pub struct ThreeAndFour;

impl ShiftPolicy for NoShift {}

impl ShiftPolicy for ThreeWayOnly {
    #[inline(always)]
    fn three_way() -> bool { true }
}

impl ShiftPolicy for FourWayOnly {
    #[inline(always)]
    fn four_way() -> bool { true }
}

impl ShiftPolicy for ThreeAndFour {
    #[inline(always)]
    fn three_way() -> bool { true }
    #[inline(always)]
    fn four_way() -> bool { true }
}

