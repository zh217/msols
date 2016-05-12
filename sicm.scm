(define win2 (frame 0.0 :pi/2 0.0 1.2))

(define ((parametric-path-action Lagrangian t0 q0 t1 q1)
         intermediate-qs)
  (let ((path (make-path t0 q0 t1 q1 intermediate-qs)))
    ;; display path
    (graphics-clear win2)
    (plot-function win2 path t0 t1 (/ (- t1 t0) 100))
    ;; compute action
    (Lagrangian-action Lagrangian path t0 t1)))


(define (find-path Lagrangian t0 q0 t1 q1 n)
  (let ((initial-qs (linear-interpolants q0 q1 n)))
    (let ((minimizing-qs (multidimensional-minimize
                          (parametric-path-action Lagrangian t0 q0 t1 q1)
                          initial-qs)))
      (make-path t0 q0 t1 q1 minimizing-qs))))

(define ((L-harmonic m k) local)
  (let ((q (coordinate local))
        (v (velocity local)))
    (- (* 1/2 m (square v)) (* 1/2 k (square q)))))

(define (((delta eta) f) q)
  (define (g epsilon t) ((f (+ q (* epsilon eta))) t))
  (lambda (t) ((ref (D g) 0) 0 t)))

(define (f q)
  (compose (literal-function 'F (-> (UP Real (UP* Real) (UP* Real)) Real))
           (Gamma q)))

(define (g q)
  (compose (literal-function 'G (-> (UP Real (UP* Real) (UP* Real)) Real))
           (Gamma q)))

(define q (literal-function 'q (-> Real (UP Real Real))))

(define eta (literal-function 'eta (-> Real (UP Real Real))))


((((delta eta) f) q) 't)

(se ((((delta eta) f) q) 't))

(((delta eta) f) q)

(f q)

((f q) 't)

((- (((delta eta) f) q) (((delta eta) f) q)) 't)

;; Leibnizian
(- ((((delta eta) (* f g)) q) 't)
   ((+ (* (((delta eta) f) q)
          (g q))
       (* (((delta eta) g) q)
          (f q))) 't))

;; distributive wrt addition
(- ((((delta eta) (+ f g)) q) 't)
   ((+ (((delta eta) f) q)
       (((delta eta) g) q)) 't))

;; vector space
(- ((((delta eta) (* 'c f)) q) 't)
   ((* 'c (((delta eta) f) q)) 't))

(define k (literal-function 'k))

;; chain rule
(- ((((delta eta) (lambda (q) (compose k (f q)))) q) 't)
   ((* (((delta eta) f) q)
       ((lambda (q)
          (compose (D k) (f q))) q))
    't))

;; commutation with partial derivation
((- (D (((delta eta) f) q))
    (((delta eta) (lambda (q) (D (f q)))) q))
 't)

(define ((L-central-polar m V) local)
  (let ((q (coordinate local))
        (qdot (velocity local)))
    (let ((r (ref q 0))
          (phi (ref q 1))
          (rdot (ref qdot 0))
          (phidot (ref qdot 1)))
      (- (* 1/2 m
            (+ (square rdot) (square (* r phidot))))
         (V r)))))

(define ((gravitational-energy G m1 m2) r)
  (- (/ (* G m1 m2) r)))

(se ((L-central-polar
 (/ (* 'm_1 'm_2) (+ 'm_1 'm_2))
 (gravitational-energy 'G 'm_1 'm_2))
 (up 't
     (up 'r 'phi)
     (up 'rdot 'phidot))))

(se
 (((Lagrange-equations
    (L-central-polar
     (/ (* 'm_1 'm_2) (+ 'm_1 'm_2))
     (gravitational-energy 'G 'm_1 'm_2)))
   (lambda (t) (up 'a (* 'n t))))
  't))

(se
 (((Lagrange-equations
    (lambda (local)
      (let ((theta (coordinate local))
            (thetadot (velocity local)))
        (+ (* 1/2 'm 'l 'l thetadot thetadot)
           (* 'm 'g 'l (cos theta))))))
   (literal-function 'theta))
  't))

(se
 (((Lagrange-equations
    (lambda (local)
      (let ((q (coordinate local))
            (v (velocity local)))
        (let ((x (ref q 0))
              (y (ref q 1)))
          (- (* 1/2 'm (square v))
             (+ (/ (square q) 2)
                (* x x y)
                (- (* 1/3 y y y))))))))
   (up (literal-function 'x)
       (literal-function 'y)))
  't))

(se
 (((Lagrange-equations
    (lambda (local)
      (let ((q (coordinate local))
            (qdot (velocity local)))
        (let ((theta (ref q 0))
              (phi (ref q 1))
              (alpha (ref qdot 0))
              (beta (ref qdot 1)))
          (* 1/2 'm 'R 'R
             (+ (* alpha alpha)
                (* beta beta (expt (sin theta) 2))))))))
   (up (literal-function 'theta)
       (literal-function 'phi)))
  't))

(define ((Lagrange-equations-prime Lagrangian order) q)
  (define (recur k)
    (if (= k 0)
        0
        (+ (recur (- k 1))
           (* (expt -1 (- k 1))
              ((expt D (- k 1))
               (compose ((partial k) Lagrangian)
                        (Gamma q (+ order 1))))))))
  (recur order))

(define ((H-Lagrangian m k) local)
  (let ((x (ref local 1))
        (v (ref local 2))
        (a (ref local 3)))
    (* (- 1/2)
       (+ (* m x a)
          (* k x x)))))
(se
 (((Lagrange-equations-prime (H-Lagrangian 'm 'k) 3)
   (literal-function 'x))
  't))

(define ((L-central-rectangular m U) local)
  (let ((q (coordinate local))
        (v (velocity local)))
    (- (* 1/2 m (square v))
       (U (sqrt (square q))))))

(se
 ((L-central-rectangular 'm (literal-function 'U))
  (up 't
      (up 'x 'y 'z)
      (up 'xdot 'ydot 'zdot))))

(define (r3->p3 local)
  (let ((q (coordinate local)))
    (let ((x (ref q 0))
          (y (ref q 1))
          (z (ref q 2)))
      (up (sqrt (square q))
          (acos (/ z (sqrt (square q))))
          (atan y x)))))

(define (p3->r3 local)
  (let ((polar-tuple (coordinate local)))
    (let ((r (ref polar-tuple 0))
          (theta (ref polar-tuple 1))
          (phi (ref polar-tuple 2)))
      (up (* r (sin theta) (cos phi))
          (* r (sin theta) (sin phi))
          (* r (cos theta))))))

(define (L-central-polar m U)
  (compose (L-central-rectangular m U) (F->C p3->r3)))

(se
 ((L-central-polar 'm (literal-function 'U))
  (up 't
      (up 'r 'theta 'phi)
      (up 'rdot 'thetadot 'phidot))))


(define ((L-helical-bead m d h g) local)
  (let ((x (coordinate local))
        (xdot (velocity local))
        (c (/ (* 2 'pi) h)))
    (let ((v (up xdot
                 (* -1 (sin (* c x)) c xdot)
                 (* (cos (* c x)) c xdot)))
          (V (* m g (sin (* c x)))))
      (- (* 1/2 m (square v)) V))))

(se
 ((L-helical-bead 'm 'd 'h 'g)
  ((Gamma (literal-function 'x))
   't)))

(se
 (((Lagrange-equations (L-helical-bead 'm 'd 'h 'g))
   (literal-function 'x))
  't))

(define ((L-triaxial-surface m a b c) local)
  (let ((q (coordinate local))
        (qdot (velocity local)))
    (let ((theta (ref q 0))
          (phi (ref q 1))
          (thetadot (ref qdot 0))
          (phidot (ref qdot 1)))
      (let (;(xprime (* (sin theta) (cos phi)))
            ;(yprime (* (sin theta) (sin phi)))
            ;(zprime (cos theta))
            (xprimedot (+ (* (cos theta) (cos phi) thetadot)
                          (* (sin theta) (sin phi) phidot -1)))
            (yprimedot (+ (* (cos theta) (sin phi) thetadot)
                          (* (sin theta) (cos phi) phidot)))
            (zprimedot (* -1 (sin theta) thetadot)))
        (let ((v (up (* a xprimedot)
                     (* b yprimedot)
                     (* c zprimedot))))
          (* 1/2 m (square v)))))))

(se
 ((L-triaxial-surface 'm 'a 'b 'c)
  ((Gamma (up (literal-function 'theta)
              (literal-function 'phi)))
   't)))

(se
 (((Lagrange-equations (L-triaxial-surface 'm 'a 'b 'c))
   (up (literal-function 'theta)
       (literal-function 'phi)))
  't))

(define ((L-two-bar-linkage m1 m2 m3 l1 l2 g) local)
  (let ((q (coordinate local))
        (v (velocity local)))
    (let ((x2 (ref q 0))
          (y2 (ref q 1))
          (t1 (ref q 2))
          (t2 (ref q 3))
          (vx2 (ref v 0))
          (vy2 (ref v 1))
          (vt1 (ref v 2))
          (vt2 (ref v 3)))
      (let ((x1 (+ x2 (* l1 (cos t1))))
            (y1 (+ y2 (* l1 (sin t1))))
            (x3 (+ x2 (* l2 (cos t2))))
            (y3 (+ y2 (* l2 (sin t2))))
            (vx1 (+ vx2 (* l1 (sin t1) vt1 -1)))
            (vy1 (+ vy2 (* l1 (cos t1) vt1)))
            (vx3 (+ vx2 (* l2 (sin t2) vt2 -1)))
            (vy3 (+ vy2 (* l2 (cos t2) vt2))))
        (let ((v1 (up vx1 vy1))
              (v2 (up vx2 vy2))
              (v3 (up vx3 vy3)))
          (+ (* 1/2 m1 (square v1))
             (* 1/2 m2 (square v2))
             (* 1/2 m3 (square v3))
             (* -1 m1 g y1)
             (* -1 m2 g y2)
             (* -1 m3 g y3)))))))

(se
 ((L-two-bar-linkage 'm_1 'm_2 'm_3 'l_1 'l_2 'g)
  ((Gamma (up (literal-function 'x_2)
              (literal-function 'y_2)
              (literal-function 'theta_1)
              (literal-function 'theta_2)))
   't)))

(se
 (((Lagrange-equations (L-two-bar-linkage 'm_1 'm_2 'm_3 'l_1 'l_2 'g))
   (up (literal-function 'x_2)
       (literal-function 'y_2)
       (literal-function 'theta_1)
       (literal-function 'theta_2)))
  't))

(define ((L-sliding-pendulum m1 m2 l g) local)
  (let ((q (Q local))
        (qdot (Qdot local)))
    (let ((x (ref q 0))
          (t (ref q 1))
          (vx (ref qdot 0))
          (vt (ref qdot 1)))
      (let ((x2 (+ x (* l (sin t))))
            (y2 (* l (cos t) -1))
            (vx2 (+ vx (* l (cos t) vt)))
            (vy2 (* l (sin t) vt)))
        (let ((V (* m2 g y2))
              (v2 (up vx2 vy2)))
          (+ (* 1/2 m1 (square vx))
             (* 1/2 m2 (square v2))
             (- V)))))))

(se
 ((L-sliding-pendulum 'm_1 'm_2 'l 'g)
  ((Gamma (up (literal-function 'x)
              (literal-function 'theta)))
   't)))

(se
 (((Lagrange-equations (L-sliding-pendulum 'm_1 'm_2 'l 'g))
   (up (literal-function 'x)
       (literal-function 'theta)))
  't))

(define ((L-dumbbell-formal m0 m1 l) local)
  (let ((q (Q local))
        (qdot (Qdot local)))
    (let ((x0 (ref q 0))
          (y0 (ref q 1))
          (x1 (ref q 2))
          (y1 (ref q 3))
          (F (ref q 4))
          (vx0 (ref qdot 0))
          (vy0 (ref qdot 1))
          (vx1 (ref qdot 2))
          (vy1 (ref qdot 3))
          (Fdot (ref qdot 4)))
      (- (* 1/2 (+ (* m0 vx0 vx0)
                   (* m0 vy0 vy0)
                   (* m1 vx1 vx1)
                   (* m1 vy1 vy1)))
         (* 1/2 (/ F l)
            (- (+ (square (- x1 x0))
                  (square (- y1 y0)))
               (square l)))))))

(se
 ((L-dumbbell-formal 'm_0 'm_1 'l)
  ((Gamma (up (literal-function 'x_0)
              (literal-function 'y_0)
              (literal-function 'x_1)
              (literal-function 'y_1)
              (literal-function 'F)))
   't)))

(se
 (((Lagrange-equations (L-dumbbell-formal 'm_0 'm_1 'l))
   (up (literal-function 'x_0)
       (literal-function 'y_0)
       (literal-function 'x_1)
       (literal-function 'y_1)
       (literal-function 'F)))
  't))

(define ((cm->formal m0 m1) local)
  (let ((tuple (coordinate local)))
    (let ((xcm (ref tuple 0))
          (ycm (ref tuple 1))
          (theta (ref tuple 2))
          (c (ref tuple 3))
          (F (ref tuple 4)))
      (up (- xcm (* c (cos theta) (/ m1 (+ m0 m1))))
          (- ycm (* c (sin theta) (/ m1 (+ m0 m1))))
          (+ xcm (* c (cos theta) (/ m0 (+ m0 m1))))
          (+ ycm (* c (sin theta) (/ m0 (+ m0 m1))))
          F))))

(define (L-dumbbell-formal-cm m0 m1 l)
  (compose (L-dumbbell-formal m0 m1 l)
           (F->C (cm->formal m0 m1))))

(se
 ((L-dumbbell-formal-cm 'm_0 'm_1 'l)
  ((Gamma (up (literal-function 'x_cm)
              (literal-function 'y_cm)
              (literal-function 'theta)
              (literal-function 'c)
              (literal-function 'F)))
   't)))

(se
 (((Lagrange-equations (L-dumbbell-formal-cm 'm_0 'm_1 'l))
   (up (literal-function 'x_cm)
       (literal-function 'y_cm)
       (literal-function 'theta)
       (literal-function 'c)
       (literal-function 'F)))
  't))

(define ((L-dumbbell-reduced m0 m1 l) local)
  (let ((q (Q local))
        (v (Qdot local)))
    (let ((x (ref q 0))
          (y (ref q 1))
          (t (ref q 2))
          (vx (ref v 0))
          (vy (ref v 1))
          (vt (ref v 2)))
      (let ((I (+ (* m0 (square (/ (* l m1) (+ m0 m1))))
                  (* m1 (square (/ (* l m0) (+ m0 m1)))))))
        (+ (* 1/2 (/ (* m0 m1) (+ m0 m1)) (square (up vx vy)))
           (* 1/2 I (square vt)))))))

(se
 (((Lagrange-equations (L-dumbbell-reduced 'm_0 'm_1 'l))
   (up (literal-function 'x_cm)
       (literal-function 'y_cm)
       (literal-function 'theta)))
  't))

(define ((L-driven-pendulum-formal m l g ys) local)
  (let ((t (time local))
        (q (Q local))
        (v (Qdot local)))
    (let ((x (ref q 0))
          (y (ref q 1))
          (F (ref q 2))
          (vx (ref v 0))
          (vy (ref v 1)))
      (+ (* 1/2 m (square (up vx vy)))
         (- (* m g y))
         (- (* 1/2 (/ F l) (- (square (up (- (ys t) y) x)) (square l))))))))

(se
 ((L-driven-pendulum-formal 'm 'l 'g (literal-function 'ys))
  ((Gamma (up (literal-function 'x)
              (literal-function 'y)
              (literal-function 'F)))
   't)))

(se
 (((Lagrange-equations (L-driven-pendulum-formal 'm 'l 'g (literal-function 'y_s)))
   (up (literal-function 'x)
       (literal-function 'y)
       (literal-function 'F)))
  't))

(define ((pendulum-p->c ys) local)
  (let ((t (time local))
        (q (Q local)))
    (let ((theta (ref q 0))
          (c (ref q 1))
          (F (ref q 2)))
      (up (* c (sin theta))
          (- (ys t)
             (* c (cos theta)))
          F))))

(define (L-driven-pendulum-formal-polar m l g ys)
  (compose (L-driven-pendulum-formal m l g ys)
           (F->C (pendulum-p->c ys))))

(se
 (((Lagrange-equations (L-driven-pendulum-formal-polar
                        'm 'l 'g (literal-function 'y_s)))
   (up (literal-function 'theta)
       (literal-function 'c)
       (literal-function 'F)))
  't))

(define ((L-pendulum-formal-c m l g) local)
  (let ((q (Q local))
        (v (Qdot local)))
    (let ((x (ref q 0))
          (y (ref q 1))
          (F (ref q 2))
          (vx (ref v 0))
          (vy (ref v 1)))
      (- (* 1/2 m (square (up vx vy)))
         (* m g y)
         (* 1/2 (/ F l) (- (square (up x y)) (square l)))))))

(se
 (((Lagrange-equations (L-pendulum-formal-c 'm 'l 'g))
   (up (literal-function 'x)
       (literal-function 'y)
       (literal-function 'F)))
  't))

(define (p->c local)
  (let ((q (Q local)))
    (let ((theta (ref q 0))
          (c (ref q 1))
          (F (ref q 2)))
      (up (* c (sin theta))
          (- (* c (cos theta)))
          F))))

(define (L-pendulum-formal m l g)
  (compose (L-pendulum-formal-c m l g)
           (F->C p->c)))

(se
 (((Lagrange-equations (L-pendulum-formal 'm 'l 'g))
   (up (literal-function 'theta)
       (literal-function 'c)
       (literal-function 'F)))
  't))

(define ((L-foucault-np-cart G M m R) local)
  (let ((q (Q local))
        (v (Qdot local)))
    (let ((x (ref q 0))
          (y (ref q 1))
          (z (ref q 2)))
      (+ (* 1/2 m (square v))
         ;(/ (* G M m z) (square R))
         (/ (* G M m) (sqrt (square q)))
         ))))

(se
 ((L-foucault-np-cart 'G 'M 'm 'R)
  (up 't
      (up 'x 'y 'z)
      (up 'xdot 'ydot 'zdot))))

(define ((Rz omega) tuple)
  (let ((x (ref tuple 0))
        (y (ref tuple 1))
        (z (ref tuple 2)))
    (up (+ (* x (cos omega)) (* y (sin omega) -1))
        (+ (* x (sin omega)) (* y (cos omega)))
        z)))

(define ((Ry omega) tuple)
  (let ((x (ref tuple 0))
        (y (ref tuple 1))
        (z (ref tuple 2)))
    (up (+ (* x (cos omega)) (* z (sin omega)))
        y
        (+ (* x (sin omega) -1) (* z (cos omega))))))

(se
 ((compose (Rz (* 'Omega 't)) (Ry 'phi))
  (up 'x_0 'y_0 'z_0)))

(define ((north-pole-polar->cartesian R l) tuple)
  (let ((theta (ref tuple 0))
        (lamb (ref tuple 1)))
    (up (* l (sin theta) (cos lamb))
        (* l (sin theta) (sin lamb))
        (+ R l (- (* l (cos theta)))))))

(define ((foucault-transform Omega phi R l) local)
  (let ((t (time local))
        (q (Q local)))
    (let ((theta (ref q 0))
          (lamb (ref q 1)))
      ((compose (Rz (* Omega t))
                (Ry phi))
       ((north-pole-polar->cartesian R l)
        (up theta lamb))))))

(se
 ((foucault-transform 'Omega 'phi 'R 'l)
  (up 't
      (up 'theta 'lambda)
      (up 'thetadot 'lambdadot))))

(define (L-foucault Omega phi R l M m G)
  (compose (L-foucault-np-cart G M m R)
           (F->C (foucault-transform Omega phi R l))))

(se
 ((L-foucault 'Omega 0 'R 'l 'M 'm 'G)
  (up 't
      (up 'theta 'lambda)
      (up 'thetadot 'lambdadot))))

(se
 (((Lagrange-equations (L-foucault 'Omega 0 'R 'l 'M 'm 'G))
   (up (literal-function 'theta)
       (literal-function 'lambda)))
  't))


(define ((T-pend m l g ys) local)
  (let ((t (time local))
        (theta (coordinate local))
        (thetadot (velocity local)))
    (let ((vys (D ys)))
      (* 1/2 m
         (+ (square (* l thetadot))
            (square (vys t))
            (* 2 l (vys t) thetadot (sin theta)))))))

(define ((V-pend m l g ys) local)
  (let ((t (time local))
        (theta (coordinate local)))
    (* m g (- (ys t) (* l (cos theta))))))

(define L-pend (- T-pend V-pend))

(define ((V-central alpha beta) r)
  (- (* beta (expt r alpha))))

(define ((L-central m alpha beta) local)
  (let ((q (Q local))
        (v (Qdot local)))
    (let ((r (ref q 0))
          (theta (ref q 1))
          (rdot (ref v 0))
          (thetadot (ref v 1)))
      (+ (* 1/2 m (+ (square rdot)
                     (square (* r thetadot))))
         ((V-central alpha beta) r)))))

(se
 ((L-central 'm 'alpha 'beta)
  (up 't
      (up 'r 'theta)
      (up 'rdot 'thetadot))))

(define (central-field-state-derivative m alpha beta)
  (Lagrangian->state-derivative
   (L-central m alpha beta)))

(se
 ((central-field-state-derivative 'm 'alpha 'beta)
  (up 't
      (up 'r 'theta)
      (up 'rdot 'thetadot))))

(define plot-win (frame -10.0 10.0 -10.0 10.0))

(define ((monitor-polar win) state)
  (let ((q (Q state)))
    (let ((r (ref q 0))
          (theta (ref q 1)))
      (plot-point win (* r (cos theta)) (* r (sin theta))))))

((evolve central-field-state-derivative
         1.0 ; m
         0.25 ; alpha
         2.0 ; beta
         )
 (up 0.0
     (up 3.0 0.0)
     (up 0.3 0.3))
 (monitor-polar plot-win)
 0.01
 100.0
 1.0e-13)


(se
 ((L-foucault 'Omega 'phi 'R 'l 'M 'm 'G)
  (up 't
      (up 'theta 'lambda)
      (up 'thetadot 'lambdadot))))

(define plot-win (frame 0.0 100.0 :-pi :pi))

(define (foucault-state-derivative Omega phi R l M m G)
  (Lagrangian->state-derivative
   (L-foucault Omega phi R l M m G)))

(se
 ((foucault-state-derivative 'Omega 'phi 'R 'l 'M 'm 'G)
  (up 't
      (up 'theta 'lambda)
      (up 'thetadot 'lambdadot))))

(define ((monitor-lambda win) state)
  (let ((lamb ((principal-value :pi) (ref (coordinate state) 1))))
    (plot-point win (time state) lamb)))

((evolve foucault-state-derivative
         0.1 ;Omega
         0.0 ;phi
         10.0 ;R
         0.1 ;l
         10 ;M
         0.1 ;m
         0.01 ;G
         )
 (up 0.0
     (up 0.1 0.0)
     (up 0.0 0.0))
 (monitor-lambda plot-win)
 0.01
 100.0
 1.0e-13)

(define ((L-uniform-acceleration m g) local)
  (let ((q (Q local))
        (v (velocity local)))
    (let ((z (ref q 2)))
      (- (* 1/2 m (square v)) (* m g z)))))

(define ((dp-coordinates l zs) local)
  (let ((t (time local))
        (q (Q local)))
    (let ((theta (ref q 0))
          (phi (ref q 1)))
      (let ((x (* l (sin theta) (cos phi)))
            (y (* l (sin theta) (sin phi)))
            (z (* l (cos theta))))
        (up x
            y
            (- (zs t) z))))))

(define ((L-pend m l g zs) local)
  (let ((t (time local))
        (q (Q local))
        (v (Qdot local)))
    (let ((theta (ref q 0))
          (thetadot (ref v 0)))
      (let ((total (- (* l m (cos theta) (((expt D 2) zs) t))
                      (* l m (sin theta) ((D zs) t) thetadot))))
        (+ ((compose (L-uniform-acceleration m g)
                     (F->C (dp-coordinates l zs)))
            local)
           total)))))

(se
 ((L-pend 'm 'l 'g (literal-function 'z_s))
  (up 't
      (up 'theta 'phi)
      (up 'thetadot 'phidot))))


(se
 ((Lagrangian->energy (L-pend 'm 'l 'g (literal-function 'z_s)))
  (up 't
      (up 'theta 'phi)
      (up 'thetadot 'phidot))))

(se
 (((partial 1) (L-pend 'm 'l 'g (literal-function 'z_s)))
  (up 't
      (up 'theta 'phi)
      (up 'thetadot 'phidot))))

(se
 (((partial 2) (L-pend 'm 'l 'g (literal-function 'z_s)))
  (up 't
      (up 'theta 'phi)
      (up 'thetadot 'phidot))))

(define ((L0 m V) local)
  (let ((t (time local))
        (q (coordinates local))
        (v (velocities local)))
    (- (* 1/2 m (square v)) (V t q))))

(define ((V a GM0 GM1 m) t xy)
  (let ((Omega (sqrt (/ (+ GM0 GM1) (expt a 3))))
        (a0 (* (/ GM1 (+ GM0 GM1)) a))
        (a1 (* (/ GM0 (+ GM0 GM1)) a)))
    (let ((x (ref xy 0)) (y (ref xy 1))
          (x0 (* -1 a0 (cos (* Omega t))))
          (y0 (* -1 a0 (sin (* Omega t))))
          (x1 (* +1 a1 (cos (* Omega t))))
          (y1 (* +1 a1 (sin (* Omega t)))))
      (let ((r0
             (sqrt (+ (square (- x x0)) (square (- y y0)))))
            (r1
             (sqrt (+ (square (- x x1)) (square (- y y1))))))
        (- (+ (/ (* GM0 m) r0) (/ (* GM1 m) r1)))))))

(se
 ((V 'a 'GM_0 'GM_1 'm)
  't
  (up 'x 'y)))

(define ((LR3B m a GM0 GM1) local)
  (let ((q (coordinates local))
        (qdot (velocities local))
        (Omega (sqrt (/ (+ GM0 GM1) (expt a 3))))
        (a0 (* (/ GM1 (+ GM0 GM1)) a))
        (a1 (* (/ GM0 (+ GM0 GM1)) a)))
    (let ((x (ref q 0))     (y (ref q 1))
          (xdot (ref qdot 0)) (ydot (ref qdot 1)))
      (let ((r0 (sqrt (+ (square (+ x a0)) (square y))))
            (r1 (sqrt (+ (square (- x a1)) (square y)))))
        (+ (* 1/2 m (square qdot))
           (* 1/2 m (square Omega) (square q))
           (* m Omega (- (* x ydot) (* xdot y)))
           (/ (* GM0 m) r0) (/ (* GM1 m) r1))))))

(define ((LR3B1 m a0 a1 Omega GM0 GM1) local)
  (let ((q (coordinates local))
        (qdot (velocities local)))
    (let ((x (ref q 0))     (y (ref q 1))
          (xdot (ref qdot 0)) (ydot (ref qdot 1)))
      (let ((r0 (sqrt (+ (square (+ x a0)) (square y))))
            (r1 (sqrt (+ (square (- x a1)) (square y)))))
        (+ (* 1/2 m (square qdot))
           (* 1/2 m (square Omega) (square q))
           (* m Omega (- (* x ydot) (* xdot y)))
           (/ (* GM0 m) r0) (/ (* GM1 m) r1))))))

(se
 ((Lagrangian->energy (LR3B1 'm 'a_0 'a_1 'Omega 'GM_0 'GM_1))
  (up 't (up 'x_r 'y_r) (up 'v_r^x 'v_r^y))))

(define ((Lr m a0 a1 Omega G M0 M1) local)
  (let ((q (Q local))
        (v (Qdot local)))
    (let ((x (ref q 0))
          (y (ref q 1))
          (vx (ref v 0))
          (vy (ref v 1)))
      (let ((r0 (sqrt (+ (square (+ x a0)) (square y))))
            (r1 (sqrt (+ (square (- x a1)) (square y)))))
        (+ (* 1/2 m (square v))
           (* 1/2 m (square Omega) (square q))
           (* m Omega (- (* x vy) (* vx y)))
           (/ (* G M0 m) r0)
           (/ (* G M1 m) r1))))))

(se
 ((Lr 'm 'a_0 'a_1 'Omega 'G 'M_0 'M_1)
  (up 't
      (up 'x_r 'y_r)
      (up 'xdot_r 'ydot_r))))

(se
 (((Lagrange-equations (Lr 'm 'a_0 'a_1 'Omega 'G 'M_0 'M_1))
   (up (literal-function 'x_r)
       (literal-function 'y_r)))
  't))

(define ((L-free m) local)
  (let ((qdot (Qdot local)))
    (* 1/2 m (square qdot))))

(define ((scale-trans a b c) local)
  (let ((q (Q local)))
    (let ((x (ref q 0))
          (y (ref q 1))
          (z (ref q 2)))
      (up (* x a)
          (* y b)
          (* z c)))))

(define (p->c local)
  (let ((q (Q local)))
    (let ((theta (ref q 0))
          (phi (ref q 1)))
      (up (* (sin theta) (cos phi))
          (* (sin theta) (sin phi))
          (cos theta)))))

(define (L-ellipsoid m a b c)
  (compose (L-free m)
           (F->C (scale-trans a b c))
           (F->C p->c)))

(se
 ((L-ellipsoid 'm 'a 'b 'c)
  (up 't
      (up 'theta 'phi)
      (up 'thetadot 'phidot))))

(define ((F-tilde a b c) s)
  (compose
   (Rx s)
   Q
   (F->C (scale-trans a b c))
   (F->C p->c)))

(se
 (((F-tilde 'a 'b 'c) 's)
  (up 't
      (up 'theta 'phi)
      (up 'thetadot 'phidot))))

(se
 (((D (F-tilde 'a 'b 'c)) 's)
  (up 't
      (up 'theta 'phi)
      (up 'thetadot 'phidot))))

(se
 ((* (ref (D (L-free 'm 'a 'a 'c)) 2)
     ((D (F-tilde 'a 'a 'c)) 0))
  (up 't
      (up 'theta 'phi)
      (up 'thetadot 'phidot))))

(se
 ((ref (D (L-ellipsoid 'm 'a 'a 'c)) 2)
  (up 't
      (up 'theta 'phi)
      (up 'thetadot 'phidot))))


(se
 ((Gamma-bar (lambda (q)
               (let ((r (ref q 0))
                     (theta (ref q 1)))
                 (up (* r (cos theta))
                     (* r (sin theta))))))
  (up 't
      (up 'r 'theta)
      (up 'rdot 'thetadot))))

(p->r
 ((Gamma (up (literal-function 'r)
             (literal-function 'theta)))
  't))

((Gamma (up (literal-function 'r)
            (literal-function 'theta))))

(define (p->rr fbar)
  )

(define ((L-free m) local)
  (let ((qdot (Qdot local)))
    (* 1/2 m (square qdot))))

(define ((constraint-triaxial a b c) tuple)
  (let ((x (ref tuple 0))
        (y (ref tuple 1))
        (z (ref tuple 2)))
    (- 1 (+ (/ (square x) (square a))
            (/ (square y) (square b))
            (/ (square z) (square c))))))

(define ((L-triaxial m a b c) local)
  (let ((t (time local))
        (q (Q local))
        (qdot (Qdot local)))
    (let ((x (ref q 0))
          (y (ref q 1))
          (z (ref q 2))
          (l (ref q 3))
          (vx (ref qdot 0))
          (vy (ref qdot 1))
          (vz (ref qdot 2)))
      (+ ((L-free m) (up t
                         (up x y z)
                         (up vx vy vz)))
         (* l ((constraint-triaxial a b c) (up x y z)))))))

(se
 ((L-triaxial 'm 'a 'b 'c)
  ((Gamma (up (literal-function 'x)
              (literal-function 'y)
              (literal-function 'z)
              (literal-function 'lambda)))
   't)))

(se
 (((Lagrange-equations (L-triaxial 'm 'a 'b 'c))
   (up (literal-function 'x)
       (literal-function 'y)
       (literal-function 'z)
       (literal-function 'lambda)))
  't))

(define ((L-golf m g h) local)
  (let* ((q (Q local))
         (x (ref q 0))
         (y (ref q 1))
         (z (ref q 2))
         (l (ref q 3))
         (v (Qdot local))
         (vx (ref v 0))
         (vy (ref v 1))
         (vz (ref v 2)))
    (+ (* 1/2 m (square (up vx vy vz)))
       (* -1 m g z)
       (* l (- z (h x y))))))

(se
 ((L-golf 'm 'g (literal-function 'h (-> (X Real Real) Real)))
  (up 't
      (up 'x 'y 'z 'lambda)
      (up 'xdot 'ydot 'zdot 'lambdadot))))

(se
 (((Lagrange-equations (L-golf 'm 'g (literal-function 'h (-> (X Real Real) Real))))
   (up (literal-function 'x)
       (literal-function 'y)
       (literal-function 'z)
       (literal-function 'lambda)))
  't))

(define (Rx t)
  (down (up 1 0 0)
        (up 0 (cos t) (sin t))
        (up 0 (- (sin t)) (cos t))))

(define (Rz t)
  (down (up (cos t) (sin t) 0)
        (up (- (sin t)) (cos t) 0)
        (up 0 0 1)))

(define (M theta phi psi)
  (* (Rz phi)
     (Rx theta)
     (Rz psi)))

(se
 (M 'theta 'phi 'psi))

(se ((literal-function 'f) 't))

(se ((M (literal-function 'theta)
        (literal-function 'phi)
        (literal-function 'psi))
     't))

(se
 ((D (M (literal-function 'theta)
        (literal-function 'phi)
        (literal-function 'psi)))
  't))

(se
 ((D (Euler->M (up (literal-function 'theta)
                   (literal-function 'phi)
                   (literal-function 'psi))))
  't))

(define diff
  (/ (- (D (Euler->M (up (literal-function 'theta)
                         (literal-function 'phi)
                         (literal-function 'psi))))
        (* (Euler->M (up (literal-function 'theta)
                         (literal-function 'phi)
                         (literal-function 'psi)))
           (A (up 'omega^a
                  'omega^b
                  'omega^c))))
     (sin (literal-function 'theta))))

(se
 ((* (Euler->M (up (literal-function 'theta)
                   (literal-function 'phi)
                   (literal-function 'psi)))
     (A (up 'omega^a
            'omega^b
            'omega^c)))
  't))

(se
 (/ (ref (diff 't) 2 2) ((sin (literal-function 'theta)) 't)))

(se (ref (diff 't) 2 1))

(se (/ (+ (ref (diff 't) 2 1)
          ((* (ref (diff 't) 2 2)
              (/ (* (cos (literal-function 'psi))
                    (cos (literal-function 'theta)))
                 (sin (literal-function 'theta))))
           't))
       ((* (sin (literal-function 'psi)))
        't)))


(define (A omega)
  (let ((x (ref omega 0))
        (y (ref omega 1))
        (z (ref omega 2)))
    (matrix-by-rows
     (list 0 (- z) y)
     (list z 0 (- x))
     (list (- y) x 0))))
