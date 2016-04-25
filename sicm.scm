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
