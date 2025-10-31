use crate::rules::RewriteRule;
use crate::rules::{r11::R11, 
                   r12::R12,
                   r21::R21,
                   r22::R22,
                   r23::R23,
                   r24::R24,
                   r25::R25,
                   r26::R26,
                   r27::R27,
                   r28::R28};
use once_cell::sync::Lazy;

static ALL_RULES: Lazy<Vec<&'static dyn RewriteRule>> = Lazy::new(|| {
    vec![
        &R11,
        &R12,
        &R21,
        &R22,
        &R23,
        &R24,
        &R25,
        &R26,
        &R27,
        &R28,
    ]
});

//static BIND: Lazy<Vec<&'static dyn RewriteRule>> = Lazy::new(|| {
//    vec![&R11]
//});
//
//static UNBIND: Lazy<Vec<&'static dyn RewriteRule>> = Lazy::new(|| {
//    vec![&R12]
//});
//
//static THREE_WAY: Lazy<Vec<&'static dyn RewriteRule>> = Lazy::new(|| {
//    vec![&R21, &R22, &R23, &R24, &R25, &R26]
//});
//
//static FOUR_WAY: Lazy<Vec<&'static dyn RewriteRule>> = Lazy::new(|| {
//    vec![&R27, &R28]
//});

static BIND_THREE_WAY: Lazy<Vec<&'static dyn RewriteRule>> = Lazy::new(|| {
    vec![&R11, &R21, &R22, &R23, &R24, &R25, &R26]
});

static BIND_THREE_WAY_FOUR_WAY: Lazy<Vec<&'static dyn RewriteRule>> = Lazy::new(|| {
    vec![&R11, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28]
});

pub fn all_rules() -> &'static [&'static dyn RewriteRule] {
    &ALL_RULES
}

pub fn bind_threeway() -> &'static [&'static dyn RewriteRule] {
    &BIND_THREE_WAY
}

pub fn bind_threeway_fourway() -> &'static [&'static dyn RewriteRule] {
    &BIND_THREE_WAY_FOUR_WAY
}


